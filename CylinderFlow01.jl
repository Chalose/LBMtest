#= 圆柱绕流(LBGK)
D2Q9模型:
7---3---6
| \ | / |
4---1---2
| / | \ |
8---5---9
障碍物边界采用(伪)修正反弹格式，四壁采用非平衡外推可消除迁移过程的四壁粒子数不守恒问题；
出现发散的主要原因:
1、来自于障碍物的修正反弹格式并不严格，为保证精度与稳定性前提下应使Nx、Ny取较大值；
2、入口出口边界条件不够良好；
3、初始时边界受回波影响较大，大Re时应减小um或使Nx、Ny取较大值减小影响。
=#
using Plots

# 环境参数(各参数均为格子单位)=========================================================================
Nx = 600                    # x向格子数
Ny = 150                    # y向格子数
obstx = Nx / 5              # 圆柱圆心坐标
obsty = Ny / 2
obstr = Ny / 10             # 圆柱半径
um = 0.1                    # 入口最大速度
Re = 1000                   # 雷诺数
ν = um * 2 * obstr / Re     # 动力粘度
ρ₀ = 1.0                    # 初始密度
tmax = 10000                # 最大求解时间

# 求解域==============================================================================================
x = LinRange(0, Nx, Nx + 1)
y = LinRange(0, Ny, Ny + 1)
X = x' .* ones(Ny + 1, Nx + 1)           # 右为x正向，下为y正向
Y = ones(Ny + 1, Nx + 1) .* y
δx = 1.0                                 # 空间步长(格子单位)
δt = 1.0                                 # 时间步长(格子单位)
M = length(0:δt:tmax)
# 障碍物蒙板
obst = (X .- obstx) .^ 2 + (Y .- obsty) .^ 2 .<= obstr^2
CI_obst = findall(obst)

# D2Q9================================================================================================
c₀ = δx / δt
cₛ = c₀ / sqrt(3)
τ = ν / (cₛ^2 * δt) + 1 / 2
ω = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
cx = c₀ * [0 1 0 -1 0 1 -1 -1 1]
cy = c₀ * [0 0 1 0 -1 1 1 -1 -1]
# 平衡态函数
function fEq!(u, v, ρ, Q)
    feq = ω[Q] * ρ .* (1.0 .+ 1 / cₛ^2 * (
        (cx[Q] * u + cy[Q] * v) +
        (cx[Q] * u + cy[Q] * v) .^ 2 * (1 / (2 * cₛ^2)) -
        (u .^ 2 + v .^ 2) * 0.5
    ))
    return feq
end

# 存储设置============================================================================================
gap = 50                                                      # 每间隔gap个时间步存储一次数据用于返回
NUM = Int(round(M / gap)) + 1
save_u = zeros(Ny + 1, Nx + 1, NUM)
save_v = zeros(Ny + 1, Nx + 1, NUM)

# 主迭代==============================================================================================
function mainFun()
    # 初始化
    ρ = ones(Ny + 1, Nx + 1) * ρ₀
    u = zeros(Ny + 1, Nx + 1)
    v = zeros(Ny + 1, Nx + 1)
    ρOut = ρ₀ * (1 - 16 * Nx * obstr * um^2 / (cₛ^2 * Ny^2 * Re))
    U = 1 / (2 * ρ₀ * ν) * cₛ^2 * (ρ₀ - ρOut) / Nx * Y .* (Ny .- Y)
    copy!(u, U)                 # 获取初始充分发展的初始速度场
    u[CI_obst] .= 0.0
    v[CI_obst] .= 0.0
    f = zeros(Ny + 1, Nx + 1, 9)
    # 分布函数初值(初始局域平衡)
    @views for Q in 1:9
        f[:, :, Q] = fEq!(u, v, ρ, Q)
    end

    # LBGK
    num = 1
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
            # 入口
            u[2:Ny, 1] = U[2:Ny, 1]
            v[2:Ny, 1] .= 0.0
            ρ[2:Ny, 1] = ρ[2:Ny, 2]
            # 出口
            u[2:Ny, Nx+1] = u[2:Ny, Nx]
            v[2:Ny, Nx+1] = v[2:Ny, Nx]
            ρ[2:Ny, Nx+1] = ρ[2:Ny, Nx]
            # 上边界
            u[Ny+1, :] .= 0.0
            v[Ny+1, :] .= 0.0
            # 下边界
            u[1, :] .= 0.0
            v[1, :] .= 0.0
        end
        # 非平衡外推
        @views for Q in [6, 2, 9]  # 入口
            f[2:Ny, 1, Q] = fEq!(u[2:Ny, 1], v[2:Ny, 1], ρ[2:Ny, 1], Q) + f[2:Ny, 2, Q] - fEq!(u[2:Ny, 2], v[2:Ny, 2], ρ[2:Ny, 2], Q)
        end
        @views for Q in [7, 4, 8]  # 出口
            f[2:Ny, Nx+1, Q] = fEq!(u[2:Ny, Nx+1], v[2:Ny, Nx+1], ρ[2:Ny, Nx+1], Q) + f[2:Ny, Nx, Q] - fEq!(u[2:Ny, Nx], v[2:Ny, Nx], ρ[2:Ny, Nx], Q)
        end
        @views for Q in [7, 3, 6]  # 下边界
            f[1, :, Q] = fEq!(u[1, :], v[1, :], ρ[1, :], Q) + f[2, :, Q] - fEq!(u[2, :], v[2, :], ρ[2, :], Q)
        end
        @views for Q in [8, 5, 9]  # 上边界
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
uv = zeros(Ny + 1, Nx + 1, NUM);
Ω = zeros(Ny + 1, Nx + 1, NUM);
u[CI_obst, :] .= 0.0;    # 刨除圆柱形状
v[CI_obst, :] .= 0.0;
# 速度模
@views for num in 1:NUM
    uv[:, :, num] = sqrt.(u[:, :, num] .^ 2 + v[:, :, num] .^ 2)
end
# 涡量场
for num in 1:NUM
    for i in 2:Ny
        for j in 2:Nx
            Ω[i, j, num] = 1 / (2 * 1) * ((v[i, j+1, num] - v[i, j-1, num]) - (u[i+1, j, num] - u[i-1, j, num]))
        end
    end
    @views begin
        Ω[:, 1, num] = Ω[:, 2, num]
        Ω[:, Nx+1, num] = Ω[:, Nx, num]
        Ω[1, :, num] = Ω[2, :, num]
        Ω[Ny+1, :, num] = Ω[Ny, :, num]
    end
end
# 动态图
begin
    anime = @animate for num in 1:NUM
        ttl = string("Re = ", Re, "\n", "t = ", (num - 1) * gap, " [Lattice unit]", "\n", "Speed  [Lattice unit]")
        p1 = Plots.heatmap(x, y, uv[:, :, num],
            color=:turbo,
            title=ttl,
            aspect_ratio=1,
            xlims=(0, Nx), ylims=(0, Ny),
            xlabel="x [Lattice unit]", ylabel="y [Lattice unit]"
        )

        p2 = Plots.heatmap(x, y, Ω[:, :, num],
            color=:berlin,
            aspect_ratio=1,
            title="Vorticity [Lattice unit]",
            xlims=(0, Nx), ylims=(0, Ny),
            clims=(-0.04, 0.04)
        )

        p3 = Plots.plot(p1, p2, layout=(2, 1), size=(900, 800))
    end
    gif(anime, fps=8)
end
