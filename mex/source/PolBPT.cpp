// 基于 Wishart 判据的 BPT 构造与剪枝
// 理论上说，当区域规模较大，或者说合并到树的高层区域时，应当使用G0分布替代Wishart分布
// 但G0分布的参数估计很麻烦，且区域规模的阈值也很难界定，因此我们不使用G0分布。
// G0_P 分布添加了 alpha
// 等参数从而可以描述完全同质->极度异质的区域。不过就广义似然比
// 的计算和实际应用而言，我们更关注散射机制，从而去判决这两个区域应不应该合并。我们不需
// 要区分这两个区域的纹理，给出它们的类别。因此，我们常认为 alpha
// 参数是一个先验的量。 在这样的情况下，G0_P 分布具有与 Wishart 分布完全相同的
// GLR 形式。
// 参见：C. C. Freitas, A. C. Frery和A. H. Correia, 《The polarimetric 𝒢
// distribution for SAR data analysis》

#include "IPSART.h"
#include "util.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <queue>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

using namespace std;
using namespace matlab;

template <typename Float, int Dim>
using Matrix = Eigen::Matrix<complex<Float>, Dim, Dim>;

using Merger = tuple<vector<int32_t>, vector<int32_t>, vector<int32_t>,
                     vector<double>, vector<double>>;

template <typename Float>
constexpr Float eps = numeric_limits<Float>::epsilon();

// 结点结构体，每个结点对应一个区域
template <typename Float, int Dim> struct Node {
  int id;                   // 结点ID
  bool active;              // 是否存活
  int N;                    // 区域的大小
  Matrix<Float, Dim> sum_C; // 区域内协方差矩阵的统计和
  double wishart_term;      // 缓存用于计算相似度指标的 Wishart 似然项
  double sum_fnorm;         // 缓存当前区域各像素协方差矩阵F范数平方和
  double phi;               // 缓存当前区域的同质化指标Phi
};

// 边结构体 (用于优先队列)
struct Edge {
  int u, v;
  double cost;

  Edge(int i, int j, double c) : u(i), v(j), cost(c) {}
  // 最小堆需要重载大于号
  bool operator>(const Edge &other) const { return cost > other.cost; }
};

// 计算对数行列式项: L = N * ln |C / N| = N * (ln |C| - 3 * ln N)
// 其中行列式的对数项可以通过 Cholesky 分解快速计算
template <typename Float, int Dim>
inline double log_det_term(const Matrix<Float, Dim> &C, int N) {
  Eigen::LLT<Matrix<Float, Dim>> llt(C);
  Matrix<Float, Dim> L = llt.matrixL();
  double ln_det = 0.0;
  for (int i = 0; i < Dim; i++) {
    ln_det += log(static_cast<double>(L(i, i).real()));
  }
  ln_det *= 2.0;
  return (ln_det - Dim * log(N)) * N;
}

// 计算两个节点合并的 GLR 代价：D = L_merge - L_i - L_j
template <typename Float, int Dim>
inline double merge_cost(const Node<Float, Dim> &n1,
                         const Node<Float, Dim> &n2) {
  int Nm = n1.N + n2.N;
  Matrix<Float, Dim> sum_Cm = n1.sum_C + n2.sum_C;
  double d = log_det_term(sum_Cm, Nm) - n1.wishart_term - n2.wishart_term;
  return d > 0.0 ? d : 0.0;
}

// BPT 树的构建
template <typename Float, int Dim>
static void bpt_build(const ipsart::PolMatView<Float, Dim> &pol_mat,
                      vector<int32_t> &us, vector<int32_t> &vs,
                      vector<int32_t> &ps, vector<double> &costs,
                      vector<double> &phis) {
  // --- 初始化 ---
  // 叶子数就是像素数
  int rows = pol_mat.rows;
  int cols = pol_mat.cols;
  int num_leaves = rows * cols;
  // BPT树是一个真二叉树，其总结点数为 2 * num_leaves - 1
  int max_nodes = 2 * num_leaves - 1;

  // 保存树的各个结点
  vector<Node<Float, Dim>> nodes(max_nodes);
  for (int i = 0; i < num_leaves; i++) {
    nodes[i].id = i;
    nodes[i].active = true;
    nodes[i].N = 1;
    nodes[i].sum_C = pol_mat.at(i);
    // 显式确保矩阵的正定性
    for (int d = 0; d < Dim; d++) {
      nodes[i].sum_C(d, d) += numeric_limits<Float>::epsilon();
    }
    // 计算似然项
    nodes[i].wishart_term = log_det_term(nodes[i].sum_C, nodes[i].N);
    // 计算F范数平方
    nodes[i].sum_fnorm = static_cast<double>(nodes[i].sum_C.squaredNorm());
    nodes[i].phi = 0.0;
  }

  // 构建图数据结构
  vector<vector<int>> adj(max_nodes); // 邻接表
  for (int i = 0; i < num_leaves; i++) {
    adj[i].reserve(4);
  }
  priority_queue<Edge, vector<Edge>, greater<Edge>> pq; // 最小堆
  for (int c = 0; c < cols; c++) {
    for (int r = 0; r < rows; r++) {
      int cur = c * rows + r;
      // 向邻接表添加元素，并将边距推入堆
      auto add_neighbor = [&](int neighbor) {
        adj[cur].push_back(neighbor);
        adj[neighbor].push_back(cur);
        double cost = merge_cost(nodes[cur], nodes[neighbor]);
        pq.emplace(cur, neighbor, cost);
      };
      // 添加 4 邻域连接
      if (c + 1 < cols) { // 右
        int right = (c + 1) * rows + r;
        add_neighbor(right);
      }
      if (r + 1 < rows) { // 下
        int down = c * rows + r + 1;
        add_neighbor(down);
      }
    }
  }

  // 设置输出数据结构
  us.resize(num_leaves - 1);
  vs.resize(num_leaves - 1);
  ps.resize(num_leaves - 1);
  costs.resize(num_leaves - 1);
  phis.resize(num_leaves - 1);

  // --- 开始合并结点 ---
  int current_node_count = num_leaves; // 每合并（循环）一次，结点数减一
  int next_id = num_leaves;            // 下一个可用 ID
  int merge_step = 0;

  // 每次迭代中所有的新邻居，放在循环外面固定占用一片内存，免得经常malloc
  vector<int> new_neighbors;
  new_neighbors.reserve(20);

  while (current_node_count > 1 && !pq.empty()) {
    // A. Pop 最小边
    Edge top = pq.top();
    pq.pop();

    int u = top.u;
    int v = top.v;

    // B. 惰性检查: 如果任一节点已不活跃，则此边无效
    if (!nodes[u].active || !nodes[v].active) {
      continue;
    }

    // C. 执行合并
    int new_id = next_id++;

    // 初始化新节点
    nodes[new_id].id = new_id;
    nodes[new_id].active = true;
    nodes[new_id].N = nodes[u].N + nodes[v].N;
    nodes[new_id].sum_C = nodes[u].sum_C + nodes[v].sum_C;
    nodes[new_id].wishart_term =
        log_det_term(nodes[new_id].sum_C, nodes[new_id].N);
    nodes[new_id].sum_fnorm = nodes[u].sum_fnorm + nodes[v].sum_fnorm;
    Matrix<Float, Dim> zx =
        nodes[new_id].sum_C / nodes[new_id].N; // 区域内平均协方差矩阵
    double n_zx_fnorm = static_cast<double>(zx.squaredNorm()) * nodes[new_id].N;
    double sum_diff_fnorm =
        nodes[new_id].sum_fnorm - n_zx_fnorm; // 一阶中心矩的F范数平方和
    nodes[new_id].phi = sum_diff_fnorm / (n_zx_fnorm + eps<double>);

    // 记录输出: 子结点u、子结点v、父节点p、合并代价cost、同质度phi
    us[merge_step] = u + 1;
    vs[merge_step] = v + 1;
    ps[merge_step] = new_id + 1;
    costs[merge_step] = top.cost;
    phis[merge_step] = nodes[new_id].phi;
    merge_step++;

    // D. 标记旧节点死亡
    nodes[u].active = false;
    nodes[v].active = false;

    // E. 更新邻居关系
    // 新邻居集合 = (Neighbors(u) U Neighbors(v)) - {u, v}
    new_neighbors.clear();
    new_neighbors.insert(new_neighbors.end(), adj[u].begin(), adj[u].end());
    new_neighbors.insert(new_neighbors.end(), adj[v].begin(), adj[v].end());
    // STL 的 Unique 函数的作用是删除相邻重复元素
    // 因此要想实现去重效果需要先排序
    sort(new_neighbors.begin(), new_neighbors.end());
    auto last = unique(new_neighbors.begin(), new_neighbors.end());
    new_neighbors.erase(last, new_neighbors.end());

    // 更新邻接表，并计算新距离、推入堆
    adj[new_id].reserve(new_neighbors.size());
    for (int nb : new_neighbors) {
      if (nb == u || nb == v || !nodes[nb].active) {
        continue;
      }
      // 更新新节点邻接表
      adj[new_id].push_back(nb);
      // 更新邻居的邻接表, u和v已死，不会再被选中
      adj[nb].push_back(new_id);

      // 计算新连接的距离，入堆
      double cost = merge_cost(nodes[new_id], nodes[nb]);
      pq.emplace(new_id, nb, cost);
    }

    // 清理旧节点的邻接表 (使用swap trick释放内存，clear函数不会释放内存)
    vector<int>().swap(adj[u]);
    vector<int>().swap(adj[v]);

    current_node_count--;

    // F. 垃圾回收
    // 当图像分辨率很高，场景均匀性较低时，会形成大量不规则的、邻居很多的小区域，从而导致pq中存在大量无效边
    // 对此，我们在pq的大小超过一定阈值时，删除其内部的死边，回收内存，提升效率
    if (pq.size() > 5ull * current_node_count) {
      vector<Edge> clean_edges;
      clean_edges.reserve(pq.size()); // 预分配一定的内存

      while (!pq.empty()) {
        Edge e = pq.top();
        pq.pop();
        // 只有当两端都存活时才保留
        if (nodes[e.u].active && nodes[e.v].active) {
          clean_edges.push_back(e);
        }
      }

      pq = priority_queue<Edge, vector<Edge>, greater<Edge>>(
          greater<Edge>(), std::move(clean_edges));
    }
  }
}

namespace ipsart {
namespace PolBPT {

vector<data::Array> Build(const vector<data::Array> &input,
                          data::ArrayFactory &af) {
  auto in_num = input.size();
  auto in_type = input[0].getType();
  auto in_dims = input[0].getDimensions();
  assert(input.size() == 4 || input.size() == 9 || input.size() == 16);
  assert(ranges::all_of(
      input | views::take(input.size() - 1), [&](const auto &arr) {
        return arr.getDimensions() == in_dims && arr.getType() == in_type;
      }));

  vector<int32_t> us, vs, ps;
  vector<double> costs, phis;

  auto execute = [&]<typename Ty>() {
    if (in_num == 4) {
      PolMatView<Ty, 2> M(input | ranges::to<vector<data::TypedArray<Ty>>>());
      bpt_build(M, us, vs, ps, costs, phis);
    } else if (in_num == 9) {
      PolMatView<Ty, 3> M(input | ranges::to<vector<data::TypedArray<Ty>>>());
      bpt_build(M, us, vs, ps, costs, phis);
    } else {
      PolMatView<Ty, 4> M(input | ranges::to<vector<data::TypedArray<Ty>>>());
      bpt_build(M, us, vs, ps, costs, phis);
    }
  };

  if (in_type == data::ArrayType::SINGLE) {
    execute.template operator()<float>();
  } else {
    execute.template operator()<double>();
  }

  vector<size_t> out_dim{1, us.size()};
  return {
      af.createArray(out_dim, us.begin(), us.end()),
      af.createArray(out_dim, vs.begin(), vs.end()),
      af.createArray(out_dim, ps.begin(), ps.end()),
      af.createArray(out_dim, costs.begin(), costs.end()),
      af.createArray(out_dim, phis.begin(), phis.end()),
  };
}

vector<data::Array> FindRoot(const vector<data::Array> &input,
                             data::ArrayFactory &af) {
  assert(input.size() == 1);

  const data::TypedArray<int32_t> map_idx = input[0];
  auto map_idx_sp = marray_to_span(map_idx);

  auto N = map_idx_sp.size();
  vector<int32_t> final_map(N);
// 路径压缩
#pragma omp parallel for
  for (int32_t i = 0; i < N; ++i) {
    int32_t curr = i;
    int32_t parent = map_idx_sp[i] - 1;
    while (parent != curr) {
      curr = parent;                 // 移动到父节点
      parent = map_idx_sp[curr] - 1; // 获取爷爷节点
    }
    final_map[i] = curr + 1;
  }

  return {af.createArray({1, N}, final_map.begin(), final_map.end())};
}

vector<data::Array> Prune(const vector<data::Array> &input,
                          data::ArrayFactory &af) {
  assert(input.size() == 4);

  const data::TypedArray<int32_t> us = input[0];
  const data::TypedArray<int32_t> vs = input[1];
  const data::TypedArray<double> phis = input[2];
  const double th = input[3][0];
  const int N = us.getNumberOfElements();

  // 返回结果
  vector<int32_t> pruned;
  pruned.reserve(N / 2);

  // 自顶向下广度搜索整棵树
  deque<int32_t> q;
  q.push_back(N);
  while (!q.empty()) {
    int32_t curr = q.front();
    q.pop_front();

    double phi = phis[curr - 1];
    // 该结点同质度大于阈值，剪除该结点
    if (phi > th) {
      pruned.push_back(curr);
      // 将它的孩子推入队列
      int32_t l = us[curr - 1];
      if (l > 0) {
        int32_t r = vs[curr - 1];
        assert(r > 0);
        q.push_back(l);
        q.push_back(r);
      }
    }
  }

  return {af.createArray({1, pruned.size()}, pruned.begin(), pruned.end())};
}

} // namespace PolBPT
} // namespace ipsart
