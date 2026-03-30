# Geometric 模块重构设计方案

## 1. 背景与目标

### 现状

`Geometric/` 目录现有 8 个文件，包含 2 个类和 6 个独立函数，职责混合在一起：

| 文件 | 类型 | 作用 |
|------|------|------|
| `OrbitVector.m` | 类 | 原始轨道数据容器（位置+速度+时间戳） |
| `OrbitPolynomial.m` | 类 | 轨道多项式拟合与求值 |
| `Times2GroundXYZ.m` | 函数 | (azimuthTime, rangeTime) → 地面 XYZ（牛顿迭代） |
| `GroundXYZ2Times.m` | 函数 | 地面 XYZ → (azimuthTime, rangeTime)（牛顿迭代逆） |
| `Pixel2RangeTime.m` | 函数 | 像素索引 → 斜距时间 |
| `RangeTime2Pixel.m` | 函数 | 斜距时间 → 像素索引 |
| `Line2AzimuthTime.m` | 函数 | 方位线号 → 方位时间 |
| `AzimuthTime2Line.m` | 函数 | 方位时间 → 方位线号 |

**存在的问题：**
- `OrbitVector` 和 `OrbitPolynomial` 职责不清，实际是一个数据容器 + 一个拟合工具，强行拆成两个类
- 所有函数平铺在同一目录，缺乏语义分组
- 轨道模型和成像几何参数散落在函数签名里，没有统一管理
- 科研场景需要独立调用，但封装层次缺失

### 重构目标

1. **分层清晰** — 按坐标空间 nature 分组到子包
2. **保留独立入口** — 所有原有独立函数接口不变，科研用户无需学习新 API
3. **Orbit 类统一管理** — 合并 `OrbitVector` + `OrbitPolynomial`，支持原始星历和预拟合多项式两种模式
4. **可选 SARGeometry helper** — 提供方便的封装，不强制使用

---

## 2. 目录结构

```
Geometric/
  +Orbit/
    Orbit.m                  <- 合并 OrbitVector + OrbitPolynomial
  +Imaging/
    Pixel2RangeTime.m        <- 独立函数
    RangeTime2Pixel.m        <- 独立函数
    Line2AzimuthTime.m        <- 独立函数
    AzimuthTime2Line.m        <- 独立函数
  +Transform/
    Times2GroundXYZ.m         <- 独立函数
    GroundXYZ2Times.m         <- 独立函数
  SARGeometry.m               <- 可选 helper 类（持有 Orbit + 成像参数）
```

**分组逻辑：**
- `Orbit/` — 卫星轨道数据与模型
- `Imaging/` — 像素/线号 ↔ 成像时间（SAR 传感器参数层面的线性映射）
- `Transform/` — 卫星姿态+斜距时间 ↔ 地面三维坐标（涉及几何迭代求解）
- `SARGeometry.m` — 可选封装层，不归入任何子包

---

## 3. Orbit 类设计

### 3.1 设计原则

- **支持两种数据模式**：`ephemeris`（原始采样）和 `polynomial`（预拟合多项式），二选一由构造函数决定
- **统一的求值接口**：`StateAt(time)` 自动识别当前模式并求值，对调用方透明
- **科研友好**：允许不拟合多项式直接使用原始数据；拟合操作显式调用 `FitPolynomial`
- **向后兼容**：`FromPolynomial` 工厂兼容现有 `OrbitPolynomial` 的所有用法

### 3.2 属性

```matlab
classdef Orbit
    properties
        % === 模式A: 原始星历数据 ===
        Time     (:,1) datetime       % 采样时刻（非相对时间，避免歧义）
        Position (:,3) double         % XYZ in ECEF [m]
        Velocity (:,3) double         % Vx Vy Vz in ECEF [m/s]

        % === 模式B: 预拟合多项式 ===
        PolynomialCoef       (:,:) double  % 行: 多项式阶数, 列: x/y/z/vx/vy/vz
        PolynomialDegree     (1,1) double
        PolynomialTimeRange  (1,2) double  % 拟合区间 [t_start, t_end]

        % === 内部状态 ===
        Ticks  (1,:) double           % 相对时间 [s]，模式A时自动计算
    end

    properties (Dependent, Abstract)
        Length                % 数据点数（模式A）
    end

    properties (Constant)
        VALID_MODES = {'ephemeris', 'polynomial'}
    end
end
```

**设计决策说明：**
- `Time` 使用 `datetime` 类型而非相对时间 `double`，避免隐式假设（以谁为基准？谁的单位是秒还是天？）
- `Position`/`Velocity` 不拆成 X/Y/Z 四数组，减少 property 数量，降低不一致风险
- `Ticks` 为 dependent property，模式A时从 `Time` 推导，模式B时来自多项式区间

### 3.3 方法

#### 工厂方法

```matlab
methods (Static)

    % 工厂1: 从原始星历数据构造
    % 输入: pos(n,3), vel(n,3), time(n,1) datetime
    % 选项: FitDegree - 构造后立即拟合多项式（默认不拟合）
    obj = FromEphemeris(pos, vel, time, options)

    % 工厂2: 从预拟合多项式系数构造（兼容现有 OrbitPolynomial 用法）
    % 输入: coef(degree+1, 6), degree, timeRange(1,2)
    obj = FromPolynomial(coef, degree, timeRange)

    % 工厂3: 从 TLE 数据构造（预留扩展接口）
    obj = FromTLE(tle_lines, varargin)
end
```

#### 实例方法

```matlab
methods

    % 统一求值接口（自动识别模式）
    state = StateAt(obj, time)
    %   模式A: 插值返回（线性插值或拉格朗日插值）
    %   模式B: 多项式求值
    %   返回 1x6 [pos vel]

    pos = PositionAt(obj, time)
    vel = VelocityAt(obj, time)

    % 多项式拟合（仅模式A可用）
    FitPolynomial(obj, degree, varargin)
    %   选项: Normalize - 是否做时间归一化（默认 true）
    %   拟合结果存入 PolynomialCoef, PolynomialDegree, PolynomialTimeRange

    % 属性访问器
    mode = get.Mode(obj)       % 返回 'ephemeris' | 'polynomial'
    len = get.Length(obj)      % 数据点数（模式A时有效）
end
```

#### 插值策略（模式A）

原始星历数据是离散采样，`StateAt` 需要插值。默认使用**线性插值**，对于高精度场景可扩展为拉格朗日插值：

```matlab
methods (Access = private)
    state = InterpolateState(obj, time)
end
```

**决策：为何不默认用拉格朗日？**
- 线性插值在轨道采样率足够高时（>1Hz）精度足够
- 拉格朗日插值在端点附近有龙格现象，科研用户需要显式选择
- 保持简单，默认行为可预测

### 3.4 使用示例

```matlab
% === 场景1: 从原始星历构造并拟合 ===
orb = Orbit.FromEphemeris(pos, vel, time);
orb.FitPolynomial(5);
pos_at_t = orb.PositionAt(t);

% === 场景2: 直接用预拟合多项式 ===
orb = Orbit.FromPolynomial(coef, 5, [t0, t1]);
state = orb.StateAt(t);

% === 场景3: 只用原始数据插值（不拟合）===
orb = Orbit.FromEphemeris(pos, vel, time);
pos_at_t = orb.PositionAt(t);  % 线性插值
```

### 3.5 向后兼容

原 `OrbitPolynomial.FromOrbitVector(orbitVector, degree)` 的用法：

```matlab
% 原用法
poly = OrbitPolynomial.FromOrbitVector(ov, 5);

% 新用法（等价）
orb = Orbit.FromEphemeris([ov.X; ov.Y; ov.Z]', [ov.Vx; ov.Vy; ov.Vz]', ov.TimeStamp);
orb.FitPolynomial(5);
```

提供**静态适配方法**保持兼容：
```matlab
methods (Static)
    % 兼容旧 OrbitVector 接口
    obj = FromOrbitVector(orbitVector, degree)
end
```

---

## 4. SARGeometry helper 类

### 4.1 设计原则

- **持有 Orbit 实例 + 成像参数**，提供便捷的坐标变换方法
- **内部委托独立函数**，不重复实现核心算法
- **可选使用**：科研用户完全可以不碰此类，直接调用独立函数
- **不自建 Orbit**：由调用方在构造时传入，支持科研灵活性

### 4.2 属性

```matlab
classdef SARGeometry
    properties
        Orbit               Orbit        % 轨道模型（必须）
        PRF                 double       % 脉冲重复频率 [Hz]
        RangeSamplingRate   double       % 距离采样率 [Hz]
        AzimuthInitTime     double       % 方位向初始时间 [s]（相对轨道起始时间）
        RangeInitTime       double       % 距离向初始时间 [s]
        CenterXYZ    (1,3) double        % 场景中心 XYZ [m]
        EcefMajorAxis       double       % WGS84 长半轴 [m]
        EcefMinorAxis       double       % WGS84 短半轴 [m]
    end

    properties (Constant)
        SPEED_OF_LIGHT = 299792458  % m/s
    end
end
```

**设计决策：椭球参数为何单独存而非用 referenceEllipsoid？**
- 原代码直接在函数签名里传 `ecefMajorAxis` / `ecefMinorAxis`，保持一致的抽象层次
- 不引入额外依赖，保持 MATLAB 兼容性

### 4.3 方法

```matlab
methods

    %% 构造
    function obj = SARGeometry(orbit, prf, rangeSR, azInitT, rngInitT, varargin)
        % varargin: centerXYZ, ecefMajorAxis, ecefMinorAxis（可选，有默认值）
        arguments
            orbit             Orbit
            prf              (1,1) double
            rangeSR          (1,1) double
            azInitT          (1,1) double
            rngInitT         (1,1) double
            varargin{:}                     % 透传给独立函数
        end
        obj.Orbit = orbit;
        obj.PRF = prf;
        obj.RangeSamplingRate = rangeSR;
        obj.AzimuthInitTime = azInitT;
        obj.RangeInitTime = rngInitT;
        % varargin 处理 centerXYZ 和椭球参数...
    end


    %% 坐标变换方法（内部调用独立函数）
    xyz = Times2GroundXYZ(obj, azimuthTime, rangeTime)
    [azT, rngT] = GroundXYZ2Times(obj, groundPos)
    line = AzimuthTime2Line(obj, azimuthTime)
    azT  = Line2AzimuthTime(obj, line)
    rngT = Pixel2RangeTime(obj, pixel)
    pixel = RangeTime2Pixel(obj, rangeTime)
    % 注意: GroundXYZ2Times 返回 [azimuthTime, rangeTime] 两个值


    %% 便捷批量变换
    % 批量像素 → 地面坐标
    XYZ = Grid2GroundXYZ(obj, lines, pixels)
end
```

**GroundXYZ2Times 返回值调整说明：**
原 `GroundXYZ2Times.m` 返回 `[azimuthTime, rangeTime]` 两个分离值，SARGeometry 方法签名与其一致。

### 4.4 使用示例

```matlab
% === 构造 SARGeometry ===
geom = SARGeometry(orbit, prf, rangeSR, azInitT, rngInitT, ...
    'CenterXYZ', centerXYZ, 'EcefMajorAxis', 6378137, 'EcefMinorAxis', 6356752);

% === 独立函数调用（完全等价，不依赖 SARGeometry）===
xyz = Times2GroundXYZ(azimuthTime, rangeTime, orbit, centerXYZ, ...
    ecefMajorAxis, ecefMinorAxis);

% === SARGeometry 封装调用 ===
xyz = geom.Times2GroundXYZ(azimuthTime, rangeTime);

% === 批量像素 → 地面坐标 ===
[lineGrid, pixelGrid] = ndgrid(1:H, 1:W);
XYZ = geom.Grid2GroundXYZ(lineGrid, pixelGrid);
```

---

## 5. 独立函数（保持不变）

### 5.1 原则

所有原有独立函数**保留原接口不变**，存放在对应子包目录下：

```
Geometric/
  +Orbit/
    Orbit.m
  +Imaging/
    Pixel2RangeTime.m    <- function rangeTime = Pixel2RangeTime(pixel, rangeInitTime, samplingFreq)
    RangeTime2Pixel.m    <- function pixel = RangeTime2Pixel(rangeTime, rangeInitTime, samplingFreq)
    Line2AzimuthTime.m    <- function azimuthTime = Line2AzimuthTime(line, azimuthInitTime, prf)
    AzimuthTime2Line.m    <- function line = AzimuthTime2Line(azimuthTime, azimuthInitTime, prf)
  +Transform/
    Times2GroundXYZ.m     <- 原接口不变
    GroundXYZ2Times.m     <- 原接口不变
  SARGeometry.m
```

**注意**：`Times2GroundXYZ` 和 `GroundXYZ2Times` 原型保持不变，但 `orbitPlyn` 参数类型从 `OrbitPolynomial` 改为 `Orbit`：

```matlab
% 原型（修改前）
function xyz = Times2GroundXYZ(azimuthTime, rangeTime, orbitPlyn, centerXYZ, ecefMajorAxis, ecefMinorAxis, options)

% 原型（修改后，参数类型更宽泛）
function xyz = Times2GroundXYZ(azimuthTime, rangeTime, orbit, centerXYZ, ecefMajorAxis, ecefMinorAxis, options)
%   orbit: Orbit 实例（支持 polynomial 模式求值）
```

`Orbit.StateAt` 的返回值格式与原 `OrbitPolynomial.GetSataliteState` 相同（1x6 [pos vel]），所以 `Times2GroundXYZ` 内部只需要把 `orbitPlyn.GetSataliteState(azimuthTime)` 替换为 `orbit.StateAt(azimuthTime)` 即可。

### 5.2 Transform/ 函数内部改动

**Times2GroundXYZ.m** — 内部替换调用：
```matlab
% 改前
sat_state = orbitPlyn.GetSataliteState(azimuthTime);

% 改后
sat_state = orbit.StateAt(azimuthTime);
```

**GroundXYZ2Times.m** — 内部替换调用：
```matlab
% 改前
sat_state = orbitPlyn.GetSataliteState(azimuthTime);

% 改后
sat_state = orbit.StateAt(azimuthTime);
```

---

## 6. 向后兼容与迁移

### 6.1 迁移路径

| 旧 API | 新 API | 兼容性 |
|--------|--------|--------|
| `OrbitVector` | `Orbit.FromEphemeris` | 静态工厂 `Orbit.FromOrbitVector` 兼容 |
| `OrbitPolynomial` | `Orbit.FromPolynomial` | 静态工厂兼容 |
| `orbitPlyn.GetSataliteState` | `orbit.StateAt` | 返回值相同，方法名不同 |
| `Geometric/Times2GroundXYZ(...)` | `Geometric/Transform/Times2GroundXYZ(...)` | 文件移动，接口不变 |

### 6.2 路径更新

`InitIPSART.m` 中 Geometric 路径需更新：
```matlab
% 改前
addpath(genpath('Geometric'));

% 改后（显式添加子包）
addpath('Geometric');
addpath('Geometric/Orbit');
addpath('Geometric/Imaging');
addpath('Geometric/Transform');
```

---

## 7. 进一步可改进方向（不纳入本次重构）

以下内容在本次重构中**不实现**，作为未来可能的扩展点记录：

1. **LatLonHeight 坐标类型** — 目前函数名本身已体现坐标类型，引入专门类型收益有限
2. **TLE 轨道解析** — `Orbit.FromTLE` 可后续添加
3. **批量网格变换缓存** — `Grid2GroundXYZ` 可引入缓存机制避免重复计算卫星状态
4. **精度选项** — 插值阶数、迭代收敛标准等可参数化
5. **地球自转校正** — 高速平台下azimuth方向的地球自转补偿（Einstein 效应）

---

## 8. 实施顺序建议

建议按以下顺序实施，降低风险：

1. **创建 `Geometric/Orbit/Orbit.m`** — 实现合并类，加上 `FromOrbitVector` 兼容工厂
2. **创建 `Geometric/Orbit/OrbitVector.m` 兼容垫片** — 指向 `Orbit`，避免立即 break
3. **更新 `Geometric/Transform/` 函数** — 改 `orbitPlyn` 参数类型，内部调用改为 `orbit.StateAt`
4. **移动文件到子包** — `Imaging/`、`Transform/`
5. **创建 `SARGeometry.m`** — 基于更新后的独立函数
6. **更新 `InitIPSART.m`** — 路径调整
7. **移除兼容垫片** — 确认无调用方后删除旧 `OrbitVector.m` 和 `OrbitPolynomial.m`

---

## 9. 验收标准

1. 现有 `Data/SARDataReader.m` 能通过兼容垫片正常工作
2. 所有独立函数接口不变，科研用户无需修改任何调用代码
3. `Orbit.FromEphemeris` → `FitPolynomial` → `StateAt` 流程测试通过
4. `Orbit.FromPolynomial` 与原 `OrbitPolynomial` 等价性测试通过
5. `SARGeometry` 方法与对应独立函数输出数值一致
