# GNSS观测分析后端214

当前仓库：https://github.com/huangjie666777-ux/gnss-scintillation-214

纯后端服务：上传一份未压缩 RINEX 3.04 观测文件和一份 SP3-c 精密星历：

- /position — 逐历元独立解算接收机 ECEF 位置与钟差（C1C 单点定位）；
- /tec — GPS 双频（C1C/C2W/L1C/L2W）电离层 TEC 监测，输出逐星逐历元
  斜向/垂直 TEC、弧编号、穿刺点地心经纬度与排除原因。
- /reflectometry — 单频（S1C，L1 C/A 信噪比）GNSS-IR 反射测高：
  由干涉振荡反演天线相位中心下方反射面高度 H，结合天线相位中心相对
  水尺零点高程 Z 输出水位 Z - H。

## 环境

- Python 3.10.12 / FastAPI 0.115.12 / NumPy 2.2.6（全部在 .venv 中）

## 运行

    .venv/bin/python -m uvicorn main:app --host 127.0.0.1 --port 8000

## 使用

    curl -s -F "rinex=@examples/obs.rnx" -F "sp3=@examples/eph.sp3" \
         http://127.0.0.1:8000/position | python3 -m json.tool

TEC 监测（站点 WGS84 经纬高 + 逐星合并码偏差 JSON，单位 ns）：

    curl -s -F "rinex=@examples/obs.rnx" -F "sp3=@examples/eph.sp3" \
         -F "station_lat_deg=30.0" -F "station_lon_deg=114.0" \
         -F "station_height_m=50.0" \
         -F "biases=$(cat examples/biases.json)" \
         http://127.0.0.1:8000/tec | python3 -m json.tool

反射测高（站点 WGS84 经纬高、天线相位中心相对水尺零点高程 Z、
反射高度搜索网格下界/上界/步长，单位米；网格最多 5001 点）：

    curl -s -F "rinex=@examples/obs_reflect.rnx" -F "sp3=@examples/eph.sp3" \
         -F "station_lat_deg=30.0" -F "station_lon_deg=114.0" \
         -F "station_height_m=50.0" -F "antenna_height_m=12.0" \
         -F "h_min_m=1.0" -F "h_max_m=8.0" -F "h_step_m=0.005" \
         http://127.0.0.1:8000/reflectometry | python3 -m json.tool

合成样例真值：H=3.6 m、Z=12.0 m，已知水位 8.4 m。

返回每个历元：状态与失败原因、ECEF 米坐标、WGS84 经纬度与椭球高、
接收机钟差（秒）、使用/排除的卫星（含排除原因）、逐星残差（米）与 RMS。

## 生成合成示例数据

    .venv/bin/python examples/make_synthetic.py   # 写 examples/obs.rnx, examples/eph.sp3
    .venv/bin/python examples/make_reflect.py     # 写 examples/obs_reflect.rnx（仅 S1C）

## 测试

    .venv/bin/python -m pytest tests -q

## 处理规则与限制

- 仅 GPS 时间系统、GPS 卫星、C1C 伪距、正常历元（flag 0）；
  接收机钟改正预应用（RCV CLOCK OFFS APPL != 0）、事件历元、非 3.04 版本、
  压缩文件等明确拒绝并指明文件与行号。
- 按头部观测类型顺序解析固定宽度字段，兼容 SYS / # / OBS TYPES 续行；
  空白或零伪距视为缺测。
- 校验日期、卫星身份、历元数量（上限 100）、非有限值与截断。
- SP3：位置 km→m，钟差 µs→s；零坐标与 999999.999999 缺失钟差被识别。
- 卫星位置用连续 8 节点 Lagrange 插值，钟差用相邻节点线性插值；
  不跨缺测节点、不外推，不可用的卫星给出排除原因。
- 每历元独立迭代最小二乘（位置 + 钟差），由伪距与卫星钟差推发射时刻，
  补偿信号飞行期间地球自转（Sagnac）；至少 4 颗有效卫星且几何满秩才求解；
  头部近似坐标只作初值，不沿用上一历元结果。
- 不足 4 星、秩亏或不收敛的历元返回原因，其余历元继续。
- 观测头允许只声明 S1C（反射测高仅需 S1C）；/position 与 /tec 各自
  检查所需观测类型，缺失时明确拒绝。全空（所有字段空白）的卫星记录
  在解析结果中保留，下游报告为缺测而不是静默消失。
- SP3 插值不跨时间节点缺口：8 节点窗口内相邻节点间隔超过名义间隔
  1.5 倍即拒绝插值并给出原因。

## TEC 解算规则

- 观测组合：码几何无关量 GF_P = (P2 - P1) - c·DCB，载波几何无关量
  GF_L = L1·λ1 - L2·λ2（周转米）；两者均等于 40.3·STEC·(1/f2² - 1/f1²)，
  载波多一个弧内常数模糊度偏移。
- 偏差符号约定：biases 给出逐星合并 P1-P2 码偏差（DCB，ns），即码观测量
  P2 相对 P1 多含 +c·DCB 硬件延迟，因此从 (P2 - P1) 中**减去** c·DCB；
  缺偏差的卫星整星标 invalid，不静默丢弃。
- 分弧：四量缺一、L1C/L2W 的 LLI 非零、相邻历元间隔 >120 s 或载波 GF
  跳变 >1 m 均断弧；各弧独立定级（不共享偏移），以弧内
  mean(GF_P - GF_L) 为偏移加到载波上；不足 3 个完整样本的弧标 invalid。
  缺测历元在输出中保留并注明原因，不静默消失。
- 斜向 TEC = GF_L(定级后) / [40.3·(1/f2² - 1/f1²)]，单位 TECU，保留负值。
- 薄壳近似：站星射线与半径 6821 km 地心球壳求交得穿刺点（输出地心
  经纬度）；垂直 TEC = 斜向 TEC × cos(射线与壳面法向夹角)。
  站点径向取为天顶，仰角 <10° 的历元标 excluded 并注明原因。
- 卫星位置复用定位链路的 8 节点 Lagrange 插值（接收历元处取值），
  不跨缺失节点、不外推；插值不可用的历元标 excluded。
- 站点坐标与偏差做有限值/范围校验；未支持的修正
  （如 RCV CLOCK OFFS APPL ≠ 0）明确拒绝。

## 精度边界

仅 C1C 伪距单点定位，无电离层/对流层/载波/多路径改正，
合成数据下亚米级；实际观测精度为米级，受大气延迟与观测噪声主导。

## 反射测高解算规则

- 模型假设：静止平面反射面（水面）、镜面反射点在天线正下方、
  不做对流层折射改正；水位 = Z - H。
- 仅用 S1C（L1 C/A 信噪比，dB-Hz）；卫星位置复用定位链路的 8 节点
  Lagrange 插值（接收历元取值，不跨缺失节点、不外推），仰角在站点
  局部 ENU 系计算，仅取 5–25 度窗口内的历元。
- 分弧：缺 S1C、无星历、相邻历元间隔 >120 s、仰角升/落反转或任何
  被过滤的历元均断弧（弧不跨过滤缺口）；弧内仰角单调。
- 弧有效性：至少 12 个点且仰角跨度至少 5 度，否则标 failed 并注明原因。
- 反演：S1C 按 10^(S1C/20) 转线性幅度；在 sin(仰角) 上拟合并减去
  二次趋势；对网格每个 H 以最小二乘拟合相位 4*pi*H*sin(仰角)/lambda1
  的正弦、余弦与常数项，取残差平方和最小的 H（同值取较小 H）；
  非等间隔采样不做 FFT。拟合秩亏或去趋势后无剩余波动的弧标 failed；
  最优值落在网格边界时给出 warning。
- 输出：逐星逐历元记录（状态与排除原因、仰角、S1C）及逐弧结果
  （起止时间、升/落方向、样本来源 S1C 与样本列表、反射高度 H、
  水位 Z-H、幅值、残差 RMS、RSS 高度曲线与警告）。

## 模块划分

- gnss_scint214/rinex.py — RINEX 3.04 解析与校验
- gnss_scint214/sp3.py — SP3-c 解析与单位换算
- gnss_scint214/interp.py — 位置 Lagrange / 钟差线性插值
- gnss_scint214/solver.py — 逐历元最小二乘定位
- gnss_scint214/tec.py — 双频几何无关组合、分弧定级、薄壳映射与穿刺点
- gnss_scint214/reflect.py — S1C 反射测高：ENU 仰角过滤、单调分弧、
  去趋势与网格搜索反演
- gnss_scint214/geodesy.py — WGS84 坐标转换与常数
- gnss_scint214/main.py — FastAPI 入口
