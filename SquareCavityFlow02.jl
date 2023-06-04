#= 对顶盖驱动方腔流场添加障碍物(LBGK)
D2Q9模型:
7---3---6
| \ | / |
4---1---2
| / | \ |
8---5---9
障碍物边界采用(伪)修正反弹格式，四壁采用非平衡外推可消除迁移过程的四壁粒子数不守恒问题；
出现发散的主要原因来自于障碍物的修正反弹格式并不严格，为保证精度与稳定性前提下应使Nx、Ny取较大值。
单位值为格子单位制。
=#
using Plots

# 环境参数(各参数均为格子单位，即δx = 1, δt = 1)===================================================
Nx = 128                    # x向格子数
Ny = 128                    # y向格子数
obstx = 3 * Nx / 4          # 圆柱圆心坐标
obsty = Ny / 2
obstr = Nx / 10             # 圆柱半径
U = 0.6                     # 顶盖速度
Re = 5000                   # 雷诺数
ν = U * Nx / Re             # 动力粘度
ρ₀ = 1.0                    # 初始密度
tmax = 5000                 # 最大求解时间

# 求解域==========================================================================================
x = LinRange(0, Nx, Nx + 1)
y = LinRange(0, Ny, Ny + 1)
X = x' .* ones(Ny + 1, Nx + 1)           # 右为x正向，下为y正向
Y = ones(Ny + 1, Nx + 1) .* y
t = 0:1:tmax
M = length(t)
# 障碍物蒙板
obst = (X .- obstx) .^ 2 + (Y .- obsty) .^ 2 .<= obstr^2
CI_obst = findall(obst)

# D2Q9============================================================================================
# c₀ = δx/δt = 1, cₛ = c₀/sqrt(3), τ = ν/(cₛ^2*δt) + 1/2
τ = 3.0 * ν + 0.5
ω = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
cx = [0 1 0 -1 0 1 -1 -1 1]
cy = [0 0 1 0 -1 1 1 -1 -1]
# 平衡态函数
function fEq!(u, v, ρ, Q)
    feq = ω[Q] * ρ .* (1 .+ 3.0 * (cx[Q] * u + cy[Q] * v) +
                       4.5 * (cx[Q] * u + cy[Q] * v) .^ 2 -
                       1.5 * (u .^ 2 + v .^ 2))
    return feq
end

# 存储设置=========================================================================================
gap = 50             # 每间隔gap个时间步存储一次数据用于返回
NUM = Int(round(M / gap)) + 1
save_u = zeros(Ny + 1, Nx + 1, NUM)
save_v = zeros(Ny + 1, Nx + 1, NUM)

# 主函数===========================================================================================
function mainFun()
    ρ = ones(Ny + 1, Nx + 1) * ρ₀
    u = zeros(Ny + 1, Nx + 1)
    v = zeros(Ny + 1, Nx + 1)
    @views u[Ny+1, :] .= U
    f = zeros(Ny + 1, Nx + 1, 9)
    # 分布函数初值(初始局域平衡)
    @views for Q in 1:9
        f[:, :, Q] = fEq!(u, v, ρ, Q)
    end

    # LBGK迭代
    num = 1          # 存储计数
    for k in 1:M
        # 包括边界粒子的碰撞步......................................................................
        @views for Q in 1:9
            f[:, :, Q] = f[:, :, Q] * (1 - 1 / τ) + 1 / τ * fEq!(u, v, ρ, Q)
        end

        # 迁移步...................................................................................
        @views begin  # 边界迁移存在粒子数不守恒，后期由非平衡外推修正
            f[:, 2:Nx+1, 2] = f[:, 1:Nx, 2]           # 左向右
            f[:, 1:Nx, 4] = f[:, 2:Nx+1, 4]           # 右向左
            f[2:Ny+1, :, 3] = f[1:Ny, :, 3]           # 下向上
            f[1:Ny, :, 5] = f[2:Ny+1, :, 5]           # 上向下
            f[2:Ny+1, 2:Nx+1, 6] = f[1:Ny, 1:Nx, 6]   # 左下向右上
            f[2:Ny+1, 1:Nx, 7] = f[1:Ny, 2:Nx+1, 7]   # 右下向左上
            f[1:Ny, 1:Nx, 8] = f[2:Ny+1, 2:Nx+1, 8]   # 右上向左下
            f[1:Ny, 2:Nx+1, 9] = f[2:Ny+1, 1:Nx, 9]   # 左上向右下
        end

        # 障碍物修正反弹格式........................................................................
        for i in CI_obst
            f[i, [2, 3, 4, 5, 6, 7, 8, 9]] = f[i, [4, 5, 2, 3, 8, 9, 6, 7]]
        end

        # 宏观量ρ, u, v刷新.........................................................................
        @views begin
            sum!(ρ, f)
            u = sum!(u, reshape(cx .* reshape(f, (Ny + 1) * (Nx + 1), 9), Ny + 1, Nx + 1, 9)) ./ ρ
            v = sum!(v, reshape(cy .* reshape(f, (Ny + 1) * (Nx + 1), 9), Ny + 1, Nx + 1, 9)) ./ ρ
        end

        # 四壁非平衡外推............................................................................
        # 确定边界宏观量
        @views begin
            u[:, [1, Nx + 1]] .= 0.0
            v[:, [1, Nx + 1]] .= 0.0
            u[1, :] .= 0.0
            v[1, :] .= 0.0
            u[Ny+1, :] .= U
            v[Ny+1, :] .= 0.0
        end
        # 非平衡外推
        @views for Q in [6, 2, 9]  # left
            f[2:Ny, 1, Q] = fEq!(u[2:Ny, 1], v[2:Ny, 1], ρ[2:Ny, 1], Q) + f[2:Ny, 2, Q] - fEq!(u[2:Ny, 2], v[2:Ny, 2], ρ[2:Ny, 2], Q)
        end
        @views for Q in [7, 4, 8]  # right
            f[2:Ny, Nx+1, Q] = fEq!(u[2:Ny, Nx+1], v[2:Ny, Nx+1], ρ[2:Ny, Nx+1], Q) + f[2:Ny, Nx, Q] - fEq!(u[2:Ny, Nx], v[2:Ny, Nx], ρ[2:Ny, Nx], Q)
        end
        @views for Q in [7, 3, 6]  # bottom
            f[1, :, Q] = fEq!(u[1, :], v[1, :], ρ[1, :], Q) + f[2, :, Q] - fEq!(u[2, :], v[2, :], ρ[2, :], Q)
        end
        @views for Q in [8, 5, 9]  # top
            f[Ny+1, :, Q] = fEq!(u[Ny+1, :], v[Ny+1, :], ρ[Ny+1, :], Q) + f[Ny, :, Q] - fEq!(u[Ny, :], v[Ny, :], ρ[Ny, :], Q)
        end

        # 保存各向速度u、v.............................................................................
        if k == (num - 1) * gap + 1
            @views save_u[:, :, num] = u
            @views save_v[:, :, num] = v
            num += 1
        end
    end

    return save_u, save_v
end
u, v = @time mainFun();

# 绘图=============================================================================================
uv = zeros(Ny + 1, Nx + 1, NUM)
# 速度模SI[m/s]
@views for num in 1:NUM
    uv[:, :, num] = sqrt.(u[:, :, num] .^ 2 + v[:, :, num] .^ 2)
end
# 动态图
begin
    anime = @animate for num in 1:NUM
        ttl = string("Re = ", Re, "\n", "t = ", (num - 1) * gap, " [Lattice unit]")
        p1 = Plots.heatmap(x, y, uv[:, :, num],
            color=:turbo,
            title=ttl,
            colorbar_title="Speed  [Lattice unit]",
            aspect_ratio=1,
            xlims=(0, Nx), ylims=(0, Ny),
            xlabel="x [Lattice unit]", ylabel="y [Lattice unit]",
            size=(600, 700))
    end
    gif(anime, fps=8)
end
