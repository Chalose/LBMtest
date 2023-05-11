#= 采用LBGK模型，演算顶盖驱动方腔流问题
速度离散模型：D2Q9:
7---3---6
| \ | / |
4---1---2
| / | \ |
8---5---9
在流场中添加静止障碍物，壁面无滑移
=#
using Plots

# 环境参数
H = 100.0      # 方腔边长SI[m]
nu = 0.0667    # 流体运动粘度SI[m²/s]
U = 0.667      # 顶盖速率SI[m/s]
ρ₀ = 1.0       # 初始流体密度SI[kg/m³]
N = 100        # 边离散节点数
δt = 1         # 时间步SI[s]
tmax = 6000    # 计算总时间SI[s]
obstX = 70     # 圆心坐标x
obstY = 50     # 圆心坐标y
obstR = 15     # 圆半径

Re = U * H / nu
println("环境雷诺数Re = ", Re)

# 求解域
x = LinRange(0, H, N)
y = LinRange(0, H, N)
#X, Y = meshgrid2(x, y)
X = x' .* ones(N)         # 右为x正向，下为y正向
Y = ones(N)' .* y
t = 1:δt:tmax
M = length(t)
# 设置固定壁蒙板
obst = (X .- obstX) .^ 2 + (Y .- obstY) .^ 2 .<= obstR^2    # 圆壁面
obst[:, [1, N]] .= true                                 # 左右壁面
obst[1, :] .= true                                      # 底壁面
CI_obst = findall(obst)                                 # 获取固定壁的笛卡尔索引序列

# LGBK:
function LGBK()
    # D2Q9:
    δx = H / (N - 1)
    c₀ = δx / δt
    cₛ = c₀ / sqrt(3)
    τ = nu / (cₛ^2 * δt) + 1 / 2
    ω = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
    cx = c₀ * [0 1 0 -1 0 1 -1 -1 1]
    cy = c₀ * [0 0 1 0 -1 1 1 -1 -1]

    # 初始化
    ρ = ones(N, N) * ρ₀
    u = zeros(N, N, 2)
    u[N, :, 1] .= U          # 顶板滑移速度
    f = zeros(N, N, 9)
    F = zeros(N, N, 9)
    feq = zeros(N, N, 9)
    # 分布函数初值(局域平衡态):
    for k in 1:9
        f[:, :, k] = ω[k] * ρ[:, :] .* (1 .+ 1 / cₛ^2 * (
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) +
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) .^ 2 * (1 / (2 * cₛ^2)) -
            (u[:, :, 1] .^ 2 + u[:, :, 2] .^ 2) * (1 / 2)
        ))
    end

    # LGBK:
    save_P = zeros(N, N, M)
    save_ux = zeros(N, N, M)
    save_uy = zeros(N, N, M)
    # 显式步进
    for tk in 1:M
        # 碰撞步
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
        # 固定壁边界条件
        F = f
        for i in CI_obst
            f[i, 2] = F[i, 4]
            f[i, 3] = F[i, 5]
            f[i, 4] = F[i, 2]
            f[i, 5] = F[i, 3]
            f[i, 6] = F[i, 8]
            f[i, 7] = F[i, 9]
            f[i, 8] = F[i, 6]
            f[i, 9] = F[i, 7]
        end
        # 移动壁边界条件
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

    return save_P, save_ux, save_uy
end
save_P, save_ux, save_uy = @time LGBK()

# 绘图
# Plots:
umod = zeros(N, N, M)
dynamic_P = save_P
for tk in 1:M
    umod[:, :, tk] = sqrt.(save_ux[:, :, tk] .^ 2 + save_uy[:, :, tk] .^ 2)  # 速度模
    for i in CI_obst
        dynamic_P[i, tk] = 0.0                                               # 刨除障碍物的动压分布
    end
end

begin
    @gif for tk in 1:20:M
        ttl = string("t = ", round(tk, digits=1), "s ", "speed[m/s]")
        p1 = heatmap(x, y, umod[:, :, tk],
            color=:turbo,
            title=ttl,
            aspect_ratio=1)
        p2 = contour(x, y, dynamic_P[:, :, tk],
            aspect_ratio=1,
            title="dynamic pressure[pa]")
        p3 = plot(p1, p2, layout=(2, 1), size=(500, 800), xlims=(0,100), ylims=(0,100))
    end
end
