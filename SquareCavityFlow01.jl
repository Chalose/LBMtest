#= 采用LBGK模型，演算顶盖驱动方腔流问题
速度离散模型：D2Q9:
7---3---6
| \ | / |
4---1---2
| / | \ |
8---5---9
程序基于Syslab编写，默认绘图库TyPlot
=#
#using Plots
using TyPlot         # Syslab同元绘图库

# 环境参数
const H = 100.0      # 方腔边长SI[m]
const nu = 0.0667    # 流体运动粘度SI[m²/s]
const U = 0.667      # 顶盖速率SI[m/s]
const ρ₀ = 1.0       # 初始流体密度SI[kg/m³]
N = 100              # 边离散节点数
δt = 1               # 时间步SI[s]
tmax = 10000         # 计算总时间SI[s]

Re = U * H / nu
println("环境雷诺数Re = ", Re)

# LGBK
function LGBK()
    # 求解域
    x = LinRange(0, H, N)
    y = LinRange(0, H, N)
    #X, Y = meshgrid2(x, y)    # 右为x正向，下为y正向(该功能由TyPlot可提供)
    X = x' .* ones(N)         # 定义右为x正向，下为y正向
    Y = ones(N)' .* y
    t = 1:δt:tmax
    M = length(t)

    # D2Q9:
    δx = H / (N - 1)          # 空间步SI[m]
    c₀ = δx / δt
    cₛ = c₀ / sqrt(3)          # 声速
    τ = nu / (cₛ^2 * δt) + 1 / 2   # 无量纲松弛时间    
    ω = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]   # 权系数
    cx = c₀ * [0 1 0 -1 0 1 -1 -1 1]                        # 离散分子速度x分量
    cy = c₀ * [0 0 1 0 -1 1 1 -1 -1]                        # 离散分子速度y分量

    # 初始化
    ρ = ones(N, N) * ρ₀      # 初始宏观密度分布
    u = zeros(N, N, 2)       # 初始宏观速度(1、2分布为x、y方向)
    u[N, :, 1] .= U          # 顶板滑移速度
    f = zeros(N, N, 9)       # 分布函数
    feq = zeros(N, N, 9)     # 平衡分布函数
    # 分布函数初值(局域平衡态):
    for k in 1:9
        f[:, :, k] = ω[k] * ρ[:, :] .* (1 .+ 1 / cₛ^2 * (
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) +
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) .^ 2 * (1 / (2 * cₛ^2)) -
            (u[:, :, 1] .^ 2 + u[:, :, 2] .^ 2) * (1 / 2)
        ))
    end

    # LGBK:
    save_P = zeros(N, N, M)     # 保存(动)压强数据
    save_ux = zeros(N, N, M)    # 保存速率数据
    save_uy = zeros(N, N, M)
    # 显式步进:
    for tk in 1:M
        # 碰撞步:
        for k in 1:9
            feq[:, :, k] = ω[k] * ρ[:, :] .* (1 .+ 1 / cₛ^2 * (
                (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) +
                (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) .^ 2 * (1 / (2 * cₛ^2)) -
                (u[:, :, 1] .^ 2 + u[:, :, 2] .^ 2) * (1 / 2)
            ))
            f[:, :, k] = f[:, :, k] * (1 - 1 / τ) + 1 / τ * feq[:, :, k]
        end
        # 迁移步:
        f[:, 2:N, 2] = f[:, 1:N-1, 2]         # 左向右
        f[:, 1:N-1, 4] = f[:, 2:N, 4]         # 右向左
        f[2:N, :, 3] = f[1:N-1, :, 3]         # 下向上
        f[1:N-1, :, 5] = f[2:N, :, 5]         # 上向下
        f[2:N, 2:N, 6] = f[1:N-1, 1:N-1, 6]   # 左下向右上
        f[2:N, 1:N-1, 7] = f[1:N-1, 2:N, 7]   # 右下向左上
        f[1:N-1, 1:N-1, 8] = f[2:N, 2:N, 8]   # 右上向左下
        f[1:N-1, 2:N, 9] = f[2:N, 1:N-1, 9]   # 左上向右下
        # 固定壁边界条件:
        f[2:N, 1, 9] = f[2:N, 1, 7]           # 左壁，向右回弹
        f[2:N, 1, 2] = f[2:N, 1, 4]
        f[2:N, 1, 6] = f[2:N, 1, 8]

        f[2:N, N, 8] = f[2:N, N, 6]           # 右壁，向左回弹
        f[2:N, N, 4] = f[2:N, N, 2]
        f[2:N, N, 7] = f[2:N, N, 9]

        f[1, :, 6] = f[1, :, 8]               # 下壁，向上回弹
        f[1, :, 3] = f[1, :, 5]
        f[1, :, 7] = f[1, :, 9]
        # 移动壁边界条件:
        ρᵗ = f[N, 2:N-1, 1] + f[N, 2:N-1, 2] + f[N, 2:N-1, 4] + 2 * (f[N, 2:N-1, 7] + f[N, 2:N-1, 3] + f[N, 2:N-1, 6])
        f[N, 2:N-1, 5] = f[N, 2:N-1, 3]       # 上壁法向回弹
        f[N, 2:N-1, 9] = f[N, 2:N-1, 7] - 2 * ω[7] / cₛ^2 * cx[7] * U * ρᵗ    # Ladd格式，碰撞动量转移
        f[N, 2:N-1, 8] = f[N, 2:N-1, 6] - 2 * ω[6] / cₛ^2 * cx[6] * U * ρᵗ    # Ladd格式，碰撞动量转移
        # 宏观量ρ 、u
        ρ = sum(f, dims=3)
        u[:, :, 1] = sum(reshape(cx .* reshape(f, N * N, 9), N, N, 9), dims=3) ./ ρ
        u[:, :, 2] = sum(reshape(cy .* reshape(f, N * N, 9), N, N, 9), dims=3) ./ ρ
        # 刷新顶盖速度
        u[N, :, 1] .= U
        u[N, :, 2] .= 0.0
        # 保存数据动压P、各向速度ux,uy
        save_P[:, :, tk] = cₛ^2 * ρ
        save_ux[:, :, tk] = u[:, :, 1]
        save_uy[:, :, tk] = u[:, :, 2]
    end

    return X, Y, save_P, save_ux, save_uy
end
X, Y, save_P, save_ux, save_uy = @time LGBK()

# 绘图
# Typlot:
begin
    t = 1:δt:tmax
    M = length(t)
    umod = zeros(N, N, M)
    for tk in 1:M
        umod[:, :, tk] = sqrt.(save_ux[:, :, tk] .^ 2 + save_uy[:, :, tk] .^ 2)
    end

    ttl = string("t = ", round(t[end], digits=1), "s    ", "speed[m/s]")

    subplot(211)
    p1 = pcolor(X, Y, umod[:, :, end])
    title(ttl)
    axis("equal", "tight")
    colormap(p1, "jet")
    colorbar(p1)

    subplot(212)
    p2 = contour(X, Y, save_P[:, :, end], 200)
    title("dynamic pressure[pa]")
    axis("equal", "tight")
    colormap(p2, "jet")
    clabel(p2, fmt="%.4f")
end
