# Geometric 模块重构实施计划

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 按设计方案重构 Geometric 模块，创建 Orbit 统一类、子包分组、SARGeometry helper 类，同时保留所有独立函数接口。

**Architecture:** 以文件移动 + 内部 API 更新为主，不改变任何函数的外部接口。Orbit 类合并原 OrbitVector + OrbitPolynomial 的功能，SARGeometry 作为可选封装层。

**Tech Stack:** MATLAB (classdef, arguments block, datetime)

---

## 实施概览

| 顺序 | 任务 | 类型 |
|------|------|------|
| 1 | 创建 Geometric/Orbit/ 子目录 | 目录 |
| 2 | 实现 Geometric/Orbit/Orbit.m | 新文件 |
| 3 | 创建 Orbit 向后兼容垫片 | 垫片 |
| 4 | 更新 Times2GroundXYZ.m（内部调用改 StateAt）| 修改 |
| 5 | 更新 GroundXYZ2Times.m（内部调用改 StateAt）| 修改 |
| 6 | 移动 Imaging 函数到 Geometric/Imaging/ | 移动 |
| 7 | 移动 Transform 函数到 Geometric/Transform/ | 移动 |
| 8 | 创建 SARGeometry.m | 新文件 |
| 9 | 更新 InitIPSART.m 路径 | 修改 |
| 10 | 验证 Data/SARDataReader 兼容性 | 验证 |

---

## Task 1: 创建子目录结构

**Files:**
- Create: `Geometric/Orbit/` (directory)
- Create: `Geometric/Imaging/` (directory)
- Create: `Geometric/Transform/` (directory)

**Step 1: 创建目录**

```bash
mkdir -p Geometric/Orbit
mkdir -p Geometric/Imaging
mkdir -p Geometric/Transform
```

**Step 2: 验证目录存在**

```bash
ls -la Geometric/
```
Expected: 输出包含 Orbit/ Imaging/ Transform/ 三个目录

---

## Task 2: 实现 Geometric/Orbit/Orbit.m

**Files:**
- Create: `Geometric/Orbit/Orbit.m`

**设计要点：**

1. **属性** — 原始星历模式与多项式模式共用同一 Orbit 实例，通过 `Mode` 属性识别
2. **Ticks** — dependent property，模式A时由 Time 推导，模式B时由 PolynomialTimeRange 推导
3. **Time 存储** — 存 datetime 类型，Ticks 存相对秒数（用于插值和拟合）
4. **Mode 检测** — 若 PolynomialCoef 非空则多项式模式，否则星历模式

```matlab
classdef Orbit

    properties
        % === 模式A: 原始星历数据 ===
        Time     (1,:) datetime       % 采样时刻
        Position (1,:) double         % XYZ (3 x n 或 n x 3)
        Velocity (1,:) double         % VxVyVz (3 x n 或 n x 3)

        % === 模式B: 预拟合多项式 ===
        PolynomialCoef       (:,:) double
        PolynomialDegree     (1,1) double = 0
        PolynomialTimeRange  (1,2) double = [0 0]
    end

    properties (Dependent)
        Length               % 数据点数（星历模式）
        Ticks                % 相对时间向量 [s]
        Mode                 % 'ephemeris' | 'polynomial'
    end

    methods

        %% 构造（私有，只能通过工厂创建）
        function obj = Orbit()
        end


        %% 实例方法: 统一求值
        function state = StateAt(obj, time)
            % 自动识别模式
            %   星历模式: 线性插值
            %   多项式模式: 多项式求值
        end

        function pos = PositionAt(obj, time)
            state = obj.StateAt(time);
            pos = state(1:3);
        end

        function vel = VelocityAt(obj, time)
            state = obj.StateAt(time);
            vel = state(4:6);
        end

        %% 实例方法: 多项式拟合（仅星历模式）
        function FitPolynomial(obj, degree)
            % M = T * A 最小二乘拟合
            % 结果存入 PolynomialCoef, PolynomialDegree, PolynomialTimeRange
        end


        %% 属性访问器
        function len = get.Length(obj)
            if strcmp(obj.Mode, 'ephemeris')
                len = size(obj.Position, 2);
            else
                len = 0;
            end
        end

        function t = get.Ticks(obj)
            if strcmp(obj.Mode, 'ephemeris')
                time_diff = obj.Time - obj.Time(1);
                t = [time_diff.Second];
            else
                t = [];
            end
        end

        function m = get.Mode(obj)
            if ~isempty(obj.PolynomialCoef)
                m = 'polynomial';
            else
                m = 'ephemeris';
            end
        end
    end


    methods (Static)

        %% 工厂1: 从原始星历构造
        function obj = FromEphemeris(pos, vel, time)
            arguments
                pos (:,:) double
                vel (:,:) double
                time (1,:) datetime
            end
            obj = Orbit();
            obj.Position = pos;
            obj.Velocity = vel;
            obj.Time = time;
        end

        %% 工厂2: 从预拟合多项式构造
        function obj = FromPolynomial(coef, degree, timeRange)
            arguments
                coef (:,:) double
                degree (1,1) double
                timeRange (1,2) double
            end
            obj = Orbit();
            obj.PolynomialCoef = coef;
            obj.PolynomialDegree = degree;
            obj.PolynomialTimeRange = timeRange;
        end

        %% 工厂3: 兼容旧 OrbitVector 用法
        function obj = FromOrbitVector(orbitVector, degree)
            arguments
                orbitVector OrbitVector
                degree (1,1) double = 5
            end
            % 将 OrbitVector 属性转为 Orbit 星历模式
            pos = [orbitVector.X; orbitVector.Y; orbitVector.Z];
            vel = [orbitVector.Vx; orbitVector.Vy; orbitVector.Vz];
            obj = Orbit.FromEphemeris(pos, vel, orbitVector.TimeStamp);
            if nargin >= 2 && degree > 0
                obj.FitPolynomial(degree);
            end
        end
    end


    methods (Access = private)

        %% 线性插值（星历模式）
        function state = InterpolateState(obj, time)
            % time: scalar datetime
            % 返回 1x6 [pos vel]
        end

        %% 多项式求值（多项式模式）
        function state = EvaluatePolynomial(obj, time)
            % time: scalar double (相对时间)
            % 返回 1x6 [pos vel]
        end
    end

end
```

**Step 2: 创建占位实现**

先写一个最小可运行的 Orbit 类框架（Mode 属性和基本结构），确保 MATLAB 能识别这个类。

**Step 3: 补充完整实现**

随后逐步实现：
1. `FromEphemeris` / `FromPolynomial` / `FromOrbitVector` 工厂
2. `get.Length` / `get.Ticks` / `get.Mode` 访问器
3. `InterpolateState`（线性插值）
4. `EvaluatePolynomial`（时间归一化 + 求值）
5. `FitPolynomial`（从星历拟合多项式）
6. `StateAt`（自动分发到插值或求值）

**Step 4: 提交**

```bash
git add Geometric/Orbit/Orbit.m
git commit -m "feat(Geometric): add unified Orbit class combining OrbitVector + OrbitPolynomial

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 3: 创建向后兼容垫片

**Files:**
- Create: `Geometric/OrbitVector.m` (兼容垫片)
- Create: `Geometric/OrbitPolynomial.m` (兼容垫片)

**Step 1: 创建 OrbitVector 兼容垫片**

`Geometric/OrbitVector.m` 改为简单包装，内部调用 Orbit：

```matlab
classdef OrbitVector < handle
    % 兼容垫片：所有属性代理到内部 Orbit 实例
    properties
        Orbit_ (Hidden) Orbit
    end

    properties
        X; Y; Z; Vx; Vy; Vz; Ticks; TimeStamp; EcefEllipsoid
    end

    methods
        function obj = OrbitVector(x, y, z, vx, vy, vz, timestamp, reference)
            if nargin == 0
                return;
            end
            pos = [x; y; z];
            vel = [vx; vy; vz];
            obj.Orbit_ = Orbit.FromEphemeris(pos, vel, timestamp);
            obj.X = x; obj.Y = y; obj.Z = z;
            obj.Vx = vx; obj.Vy = vy; obj.Vz = vz;
            obj.Ticks = obj.Orbit_.Ticks;
            obj.TimeStamp = timestamp;
            obj.EcefEllipsoid = reference;
        end
    end
end
```

**Step 2: 创建 OrbitPolynomial 兼容垫片**

`Geometric/OrbitPolynomial.m` 改为代理到 Orbit：

```matlab
classdef OrbitPolynomial < handle
    properties
        Orbit_ (Hidden) Orbit
        Coeff
        Degree
        TimeRange
    end

    methods
        function obj = OrbitPolynomial(coefficent, degree, timeRange)
            if nargin == 0
                return;
            end
            obj.Orbit_ = Orbit.FromPolynomial(coefficent, degree, timeRange);
            obj.Coeff = coefficent;
            obj.Degree = degree;
            obj.TimeRange = timeRange;
        end

        function state = GetSataliteState(obj, time, indices)
            state = obj.Orbit_.StateAt(time);
            if nargin >= 3
                state = state(indices);
            end
        end

        function T = TimeMatrix(obj, time, degree)
            arguments
                obj OrbitPolynomial
                time (:,1) double
                degree {mustBeInteger} = obj.Degree
            end
            t_norm = 4 * (time - obj.TimeRange(1)) / diff(obj.TimeRange) - 2;
            T = ones(length(t_norm), degree + 1);
            for n = 1:degree
                T(:,n + 1) = t_norm.^n;
            end
        end

        function xyz = GetSatalitePosition(obj, relativeTime)
            xyz = obj.GetSataliteState(relativeTime, [1 2 3]);
        end

        function vxyz = GetSataliteVelocity(obj, relativeTime)
            vxyz = obj.GetSataliteState(relativeTime, [4 5 6]);
        end
    end

    methods (Static)
        function obj = FromOrbitVector(orbitVector, degree)
            arguments
                orbitVector OrbitVector
                degree {mustBeInteger} = 5
            end
            obj = OrbitPolynomial(orbitVector.Coeff, orbitVector.Degree, orbitVector.TimeRange);
        end
    end
end
```

**Step 3: 提交**

```bash
git add Geometric/OrbitVector.m Geometric/OrbitPolynomial.m
git commit -m "feat(Geometric): add backward-compat shims for OrbitVector and OrbitPolynomial

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 4: 更新 Times2GroundXYZ.m

**Files:**
- Modify: `Geometric/Times2GroundXYZ.m` (内部一行修改)
- Later Move: → `Geometric/Transform/Times2GroundXYZ.m`

**Step 1: 确认改动位置**

当前代码第 21 行：
```matlab
sat_state = orbitPlyn.GetSataliteState(azimuthTime);
```

改为：
```matlab
sat_state = orbit.StateAt(azimuthTime);
```

**注意**：函数参数名 `orbitPlyn` 保持不变（对外接口不变），内部调用改为 `orbit.StateAt`。

**Step 2: 执行修改**

```matlab
% Times2GroundXYZ.m 第21行
old = 'sat_state = orbitPlyn.GetSataliteState(azimuthTime);';
new = 'sat_state = orbit.StateAt(azimuthTime);';
```

**Step 3: 提交**

```bash
git add Geometric/Times2GroundXYZ.m
git commit -m "refactor(Geometric): Times2GroundXYZ calls orbit.StateAt instead of orbitPlyn.GetSataliteState

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 5: 更新 GroundXYZ2Times.m

**Files:**
- Modify: `Geometric/GroundXYZ2Times.m` (内部一行修改)
- Later Move: → `Geometric/Transform/GroundXYZ2Times.m`

**Step 1: 确认改动位置**

当前代码第 13 行：
```matlab
sat_state = orbitPlyn.GetSataliteState(azimuthTime);
```

**Step 2: 执行修改**

同上，`orbitPlyn.GetSataliteState(azimuthTime)` → `orbit.StateAt(azimuthTime)`

**Step 3: 提交**

```bash
git add Geometric/GroundXYZ2Times.m
git commit -m "refactor(Geometric): GroundXYZ2Times calls orbit.StateAt instead of orbitPlyn.GetSataliteState

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 6: 移动 Imaging 函数到子包

**Files:**
- Move: `Geometric/Pixel2RangeTime.m` → `Geometric/Imaging/Pixel2RangeTime.m`
- Move: `Geometric/RangeTime2Pixel.m` → `Geometric/Imaging/RangeTime2Pixel.m`
- Move: `Geometric/Line2AzimuthTime.m` → `Geometric/Imaging/Line2AzimuthTime.m`
- Move: `Geometric/AzimuthTime2Line.m` → `Geometric/Imaging/AzimuthTime2Line.m`

**Step 1: 移动文件**

```bash
mv Geometric/Pixel2RangeTime.m Geometric/Imaging/
mv Geometric/RangeTime2Pixel.m Geometric/Imaging/
mv Geometric/Line2AzimuthTime.m Geometric/Imaging/
mv Geometric/AzimuthTime2Line.m Geometric/Imaging/
```

**Step 2: 验证文件存在**

```bash
ls Geometric/Imaging/
```
Expected: AzimuthTime2Line.m  Line2AzimuthTime.m  Pixel2RangeTime.m  RangeTime2Pixel.m

**Step 3: 提交**

```bash
git add Geometric/Imaging/
git mv Geometric/Pixel2RangeTime.m Geometric/Imaging/
git mv Geometric/RangeTime2Pixel.m Geometric/Imaging/
git mv Geometric/Line2AzimuthTime.m Geometric/Imaging/
git mv Geometric/AzimuthTime2Line.m Geometric/Imaging/
git commit -m "refactor(Geometric): move imaging functions to Geometric/Imaging/ subpackage

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 7: 移动 Transform 函数到子包

**Files:**
- Move: `Geometric/Times2GroundXYZ.m` → `Geometric/Transform/Times2GroundXYZ.m`
- Move: `Geometric/GroundXYZ2Times.m` → `Geometric/Transform/GroundXYZ2Times.m`

**Step 1: 移动文件**

```bash
git mv Geometric/Times2GroundXYZ.m Geometric/Transform/
git mv Geometric/GroundXYZ2Times.m Geometric/Transform/
```

**Step 2: 提交**

```bash
git add Geometric/Transform/
git commit -m "refactor(Geometric): move transform functions to Geometric/Transform/ subpackage

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 8: 创建 SARGeometry.m

**Files:**
- Create: `Geometric/SARGeometry.m`

**实现要点：**

1. 持有 Orbit 实例 + 所有成像参数
2. 构造时传入 Orbit，不自行创建
3. 方法内部**委托调用独立函数**（复用已移动到 Transform/ 的函数）
4. options 不在属性中存储，在各方法中透传

```matlab
classdef SARGeometry

    properties
        Orbit             Orbit
        PRF               double
        RangeSamplingRate double
        AzimuthInitTime   double
        RangeInitTime     double
        CenterXYZ         (1,3) double
        EcefMajorAxis     double
        EcefMinorAxis     double
    end

    properties (Constant)
        SPEED_OF_LIGHT = 299792458
    end

    methods

        function obj = SARGeometry(orbit, prf, rangeSR, azInitT, rngInitT, varargin)
            arguments
                orbit Orbit
                prf (1,1) double
                rangeSR (1,1) double
                azInitT (1,1) double
                rngInitT (1,1) double
                varargin.CenterXYZ = [0 0 0]
                varargin.EcefMajorAxis = 6378137
                varargin.EcefMinorAxis = 6356752
            end

            obj.Orbit = orbit;
            obj.PRF = prf;
            obj.RangeSamplingRate = rangeSR;
            obj.AzimuthInitTime = azInitT;
            obj.RangeInitTime = rngInitT;
            obj.CenterXYZ = varargin.CenterXYZ;
            obj.EcefMajorAxis = varargin.EcefMajorAxis;
            obj.EcefMinorAxis = varargin.EcefMinorAxis;
        end


        %% 成像参数变换
        function rngT = Pixel2RangeTime(obj, pixel)
            rngT = obj.RangeInitTime + (pixel - 1) / (2 * obj.RangeSamplingRate);
        end

        function pixel = RangeTime2Pixel(obj, rangeTime)
            pixel = 2 * (rangeTime - obj.RangeInitTime) * obj.RangeSamplingRate + 1;
        end

        function azT = Line2AzimuthTime(obj, line)
            azT = (line - 1) / obj.PRF + obj.AzimuthInitTime;
        end

        function line = AzimuthTime2Line(obj, azimuthTime)
            line = (azimuthTime - obj.AzimuthInitTime) * obj.PRF + 1;
        end


        %% 几何变换（委托 Transform/ 独立函数）
        function xyz = Times2GroundXYZ(obj, azimuthTime, rangeTime, varargin)
            % varargin: MAX_ITER, TOLERANCE
            arguments
                obj SARGeometry
                azimuthTime (1,1) double
                rangeTime (1,1) double
                varargin.MAX_ITER = 20
                varargin.TOLERANCE = 1e-16
            end
            options = namedargs2cell(varargin);
            xyz = Times2GroundXYZ(azimuthTime, rangeTime, obj.Orbit, ...
                obj.CenterXYZ, obj.EcefMajorAxis, obj.EcefMinorAxis, options{:});
        end

        function [azT, rngT] = GroundXYZ2Times(obj, groundPos, varargin)
            arguments
                obj SARGeometry
                groundPos (1,3) double
                varargin.MAX_ITER = 20
                varargin.TOLERANCE = 1e-16
            end
            options = namedargs2cell(varargin);
            [azT, rngT] = GroundXYZ2Times(groundPos, obj.Orbit, options{:});
        end


        %% 批量变换
        function XYZ = Grid2GroundXYZ(obj, lines, pixels)
            % lines, pixels: 同维度矩阵（通常由 ndgrid 产生）
            % 返回同维度 3D XYZ 矩阵
            [azT_grid, rngT_grid] = obj.LinesPixels2Times(lines, pixels);
            XYZ = obj.Times2GroundXYZ_grid(azT_grid, rngT_grid);
        end

    end


    methods (Access = private)

        function [azT_grid, rngT_grid] = LinesPixels2Times(obj, lines, pixels)
            azT_grid = obj.Line2AzimuthTime(lines);
            rngT_grid = obj.Pixel2RangeTime(pixels);
        end

        function XYZ = Times2GroundXYZ_grid(obj, azT_grid, rngT_grid)
            [m, n] = size(azT_grid);
            XYZ = zeros(m, n, 3);
            for i = 1:m
                for j = 1:n
                    XYZ(i,j,:) = obj.Times2GroundXYZ(azT_grid(i,j), rngT_grid(i,j));
                end
            end
        end
    end

end
```

**注意**：`Grid2GroundXYZ` 目前用双层循环实现，这是最简单的版本。性能优化（向量化、缓存）作为未来改进点，不在本次范围内。

**Step 2: 提交**

```bash
git add Geometric/SARGeometry.m
git commit -m "feat(Geometric): add SARGeometry helper class as optional convenience wrapper

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 9: 更新 InitIPSART.m 路径

**Files:**
- Modify: `InitIPSART.m`

**Step 1: 读取当前 InitIPSART.m 中 Geometric 相关行**

```bash
grep -n "Geometric" InitIPSART.m
```

**Step 2: 找到当前 Geometric 路径添加方式**

预期是 `addpath(genpath('Geometric'))` 或类似。

**Step 3: 修改为显式子包路径**

```matlab
% 改前
addpath(genpath('Geometric'));

% 改后
addpath('Geometric');
addpath('Geometric/Orbit');
addpath('Geometric/Imaging');
addpath('Geometric/Transform');
```

**注意**：原 `Geometric/` 下的垫片文件（`OrbitVector.m`、`OrbitPolynomial.m`）仍在 Geometric/ 根目录，MATLAB 的 addpath 可以找到它们，无需单独添加。

**Step 4: 提交**

```bash
git add InitIPSART.m
git commit -m "chore: update InitIPSART.m to add Geometric subpackages explicitly

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## Task 10: 验证兼容性

**Files:**
- Test: `Data/SARDataReader.m`

**Step 1: 读取 Data/SARDataReader.m 中 Orbit 属性**

确认 SARDataReader.Orbit 属性类型仍为 OrbitVector（通过兼容垫片可用）：

```matlab
% 在 MATLAB 中验证
run InitIPSART.m
sdr = SARDataReader(...);
orb = sdr.Orbit;          % 应返回 OrbitVector 垫片对象
pos = orb.PositionAt(t);  % 验证 OrbitVector 仍可用
```

**Step 2: 验证 Orbit.FromOrbitVector 流程**

```matlab
ov = OrbitVector(x, y, z, vx, vy, vz, timestamp, referenceEllipsoid);
orb = Orbit.FromOrbitVector(ov, 5);
state = orb.StateAt(t);    % 应返回 1x6 向量
```

**Step 3: 验证 Times2GroundXYZ 独立函数（通过垫片）**

```matlab
% 使用 OrbitPolynomial 垫片
poly = OrbitPolynomial(coef, degree, timeRange);
xyz = Times2GroundXYZ(azimuthTime, rangeTime, poly, centerXYZ, a, b);
% poly 的 GetSataliteState 方法代理到 Orbit_.StateAt
```

**Step 4: 提交最终状态**

```bash
git status
git commit -m "test(Geometric): verify backward compatibility with SARDataReader

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>"
```

---

## 验证清单

完成后确认以下均成立：

- [ ] `Geometric/Orbit/Orbit.m` 存在且可实例化
- [ ] `Orbit.FromEphemeris` → `FitPolynomial` → `StateAt` 全流程正常
- [ ] `Orbit.FromPolynomial` 与原 OrbitPolynomial 等价
- [ ] `Orbit.FromOrbitVector` 与原 OrbitVector 兼容
- [ ] `Geometric/Imaging/` 下 4 个函数可独立调用
- [ ] `Geometric/Transform/Times2GroundXYZ` 和 `GroundXYZ2Times` 可独立调用
- [ ] `Geometric/SARGeometry` 可正常构造并调用各方法
- [ ] `Data/SARDataReader` 的 Orbit 属性仍可用
- [ ] `InitIPSART.m` 执行无错误
