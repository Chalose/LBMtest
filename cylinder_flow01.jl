#= 采用LBGK模型，演算圆柱绕流问题
速度离散模型：D2Q9:
7---3---6
| \ | / |
4---1---2
| / | \ |
8---5---9
=#
using Plots

# 环境参数============================================================================================
H = 100                                    # 流域高SI[m]
W = 500                                    # 流域宽SI[m]
obstX = W / 5                              # 圆柱位置与半径
obstY = H / 2
obstR = H / 10
u̅ᵢₙ = 1.0                                  # 入口平均速度SI[m/s]
ρ₀ = 1.0                                   # 初始密度SI[kg/m³]
Re = 200                                   # 雷诺数
δx = 1                                     # 空间步长SI[m]
δt = 0.1                                   # 时间步长SI[s]
tmax = 1000                                # 计算时间SI[s]
ν = 2*obstR*u̅ᵢₙ/Re                         # 运动粘度SI[m²/s]

# 求解域==============================================================================================
Nx = Int(W / δx)                           # x向份数
Ny = Int(H / δx)                           # y向份数
x = LinRange(0, W, Nx + 1)
y = LinRange(0, H, Ny + 1)
X = x' .* ones(Ny + 1, Nx + 1)             # 右为x正向，下为y正向
Y = ones(Ny + 1, Nx + 1) .* y
t = 0:δt:tmax
M = length(t)
# 壁面蒙板
obst = (X .- obstX) .^ 2 + (Y .- obstY) .^ 2 .<= obstR^2         # 圆柱
obst[[1, Ny + 1], :] .= true                                     # 上下壁面
CI_obst = findall(obst)

# D2Q9================================================================================================
c₀ = δx / δt
cₛ = c₀ / sqrt(3)
τ = ν / (cₛ^2 * δt) + 1 / 2
ω = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
cx = c₀ * [0 1 0 -1 0 1 -1 -1 1]
cy = c₀ * [0 0 1 0 -1 1 1 -1 -1]
# 平衡分布函数
function fEq!(rho, u, v, Q)
    F = ω[Q] * rho .* (1 .+ 1 / cₛ^2 * (
        (cx[Q] * u + cy[Q] * v) +
        (cx[Q] * u + cy[Q] * v) .^ 2 * (1 / (2 * cₛ^2)) -
        (u .^ 2 + v .^ 2) * (1 / 2)
    ))
    return F
end

# 存储设置============================================================================================
gap = 50                                                      # 每间隔gap个时间步存储一次数据用于返回
NUM = Int(round(M / gap)) + 1
save_u = zeros(Ny + 1, Nx + 1, NUM)
save_v = zeros(Ny + 1, Nx + 1, NUM)

# 主迭代==============================================================================================
function mainFun()
    # 初始化:
    ρ = ones(Ny + 1, Nx + 1) * ρ₀
    u = zeros(Ny + 1, Nx + 1)
    v = zeros(Ny + 1, Nx + 1)
    U = 6 * u̅ᵢₙ / H^2 * Y[2:Ny, 1] .* (H .- Y[2:Ny, 1])        # 入口速度抛物型分布
    @views u[2:Ny, 1, 1] = U
    f = zeros(Ny + 1, Nx + 1, 9)
    # 分布函数初值
    for Q in 1:9
        f[:, :, Q] = ω[Q] * ρ .* (1 .+ 1 / cₛ^2 * (
            (cx[Q] * u + cy[Q] * v) +
            (cx[Q] * u + cy[Q] * v) .^ 2 * (1 / (2 * cₛ^2)) -
            (u .^ 2 + v .^ 2) * (1 / 2)
        ))
        #f[:, :, Q] = fEq!(ρ, u, v, Q)
    end

    # LBGK:
    num = 1        # 存储计数
    for k in 1:M
        # 确定入出口分布f...........................................................................................................................
        # 入出口宏观边界设置(对收敛影响较大，待定)
        ρ[:, Nx+1] = ρ[:, Nx]
        ρ[:, 1] = ρ[:, 2]
        u[2:Ny, 1] = U
        v[2:Ny, 1] .= 0.0
        u[:, Nx+1] = u[:, Nx]
        # 非平衡外推(Guo格式)
        @views for Q in [2, 6, 9]
            f[2:Ny, 1, Q] = fEq!(ρ[2:Ny, 1], u[2:Ny, 1], v[2:Ny, 1], Q) + f[2:Ny, 2, Q] - fEq!(ρ[2:Ny, 2], u[2:Ny, 2], v[2:Ny, 2], Q)
        end
        @views for Q in [4, 8, 7]
            f[2:Ny, Nx+1, Q] = fEq!(ρ[2:Ny, Nx+1], u[2:Ny, Nx+1], v[2:Ny, Nx+1], Q) + f[2:Ny, Nx, Q] - fEq!(ρ[2:Ny, Nx], u[2:Ny, Nx], v[2:Ny, Nx], Q)
        end

        # 内点碰撞与迁移............................................................................................................................
        @views begin
            f[2:Ny, 2:Nx, 1] = f[2:Ny, 2:Nx, 1] + 1 / τ * (fEq!(ρ[2:Ny, 2:Nx], u[2:Ny, 2:Nx], v[2:Ny, 2:Nx], 1) - f[2:Ny, 2:Nx, 1])
            f[2:Ny, 2:Nx, 2] = f[2:Ny, 1:Nx-1, 2] + 1 / τ * (fEq!(ρ[2:Ny, 1:Nx-1], u[2:Ny, 1:Nx-1], v[2:Ny, 1:Nx-1], 2) - f[2:Ny, 1:Nx-1, 2])
            f[2:Ny, 2:Nx, 3] = f[1:Ny-1, 2:Nx, 3] + 1 / τ * (fEq!(ρ[1:Ny-1, 2:Nx], u[1:Ny-1, 2:Nx], v[1:Ny-1, 2:Nx], 3) - f[1:Ny-1, 2:Nx, 3])
            f[2:Ny, 2:Nx, 4] = f[2:Ny, 3:Nx+1, 4] + 1 / τ * (fEq!(ρ[2:Ny, 3:Nx+1], u[2:Ny, 3:Nx+1], v[2:Ny, 3:Nx+1], 4) - f[2:Ny, 3:Nx+1, 4])
            f[2:Ny, 2:Nx, 5] = f[3:Ny+1, 2:Nx, 5] + 1 / τ * (fEq!(ρ[3:Ny+1, 2:Nx], u[3:Ny+1, 2:Nx], v[3:Ny+1, 2:Nx], 5) - f[3:Ny+1, 2:Nx, 5])
            f[2:Ny, 2:Nx, 6] = f[1:Ny-1, 1:Nx-1, 6] + 1 / τ * (fEq!(ρ[1:Ny-1, 1:Nx-1], u[1:Ny-1, 1:Nx-1], v[1:Ny-1, 1:Nx-1], 6) - f[1:Ny-1, 1:Nx-1, 6])
            f[2:Ny, 2:Nx, 7] = f[1:Ny-1, 3:Nx+1, 7] + 1 / τ * (fEq!(ρ[1:Ny-1, 3:Nx+1], u[1:Ny-1, 3:Nx+1], v[1:Ny-1, 3:Nx+1], 7) - f[1:Ny-1, 3:Nx+1, 7])
            f[2:Ny, 2:Nx, 8] = f[3:Ny+1, 3:Nx+1, 8] + 1 / τ * (fEq!(ρ[3:Ny+1, 3:Nx+1], u[3:Ny+1, 3:Nx+1], v[3:Ny+1, 3:Nx+1], 8) - f[3:Ny+1, 3:Nx+1, 8])
            f[2:Ny, 2:Nx, 9] = f[3:Ny+1, 1:Nx-1, 9] + 1 / τ * (fEq!(ρ[3:Ny+1, 1:Nx-1], u[3:Ny+1, 1:Nx-1], v[3:Ny+1, 1:Nx-1], 9) - f[3:Ny+1, 1:Nx-1, 9])
        end

        # 固定壁反弹边界条件.........................................................................................................................
        for i in CI_obst
            f[i, [2, 4]] = f[i, [4, 2]]
            f[i, [3, 5]] = f[i, [5, 3]]
            f[i, [6, 8]] = f[i, [8, 6]]
            f[i, [7, 9]] = f[i, [9, 7]]
        end

        # 宏观量ρ 、u、v刷新............................................................................................................................
        sum!(ρ, f)
        u = sum!(u, reshape(cx .* reshape(f, (Ny + 1) * (Nx + 1), 9), Ny + 1, Nx + 1, 9)) ./ ρ
        v = sum!(v, reshape(cy .* reshape(f, (Ny + 1) * (Nx + 1), 9), Ny + 1, Nx + 1, 9)) ./ ρ

        # 保存各向速度u、v...........................................................................................................................
        if k == (num - 1) * gap + 1
            @views save_u[:, :, num] = u
            @views save_v[:, :, num] = v
            num += 1
        end
    end

    return save_u, save_v
end
u, v = @time mainFun()

# 绘图================================================================================================
uv = zeros(Ny + 1, Nx + 1, NUM)
Ω = zeros(Ny + 1, Nx + 1, NUM)
# 速度模SI[m/s]
@views for num in 1:NUM
    uv[:, :, num] = sqrt.(u[:, :, num] .^ 2 + v[:, :, num] .^ 2)
end
# 涡量场SI[1/s]
for num in 1:NUM
    for i in 2:Ny
        for j in 2:Nx
            Ω[i, j, num] = 1 / (2 * δx) * ((v[i, j+1, num] - v[i, j-1, num]) - (u[i+1, j, num] - u[i-1, j, num]))
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
        ttl = string("Re = ", Re, "\n", "t = ", (num - 1) * gap * δt, "s", "\n", "Speed")
        p1 = Plots.heatmap(x, y, uv[:, :, num],
            color=:turbo,
            title=ttl,
            aspect_ratio=1,
            xlims=(0, W), ylims=(0, H),
            colorbar_title="m/s",
            xlabel="x (m)", ylabel="y (m)")


        p2 = Plots.heatmap(x, y, Ω[:, :, num],
            color=:berlin,
            aspect_ratio=1,
            title="Vorticity field",
            xlims=(0, W), ylims=(0, H),
            clims=(-0.3, 0.3),
            colorbar_title="1/s",
            xlabel="x (m)", ylabel="y (m)")

        p3 = Plots.plot(p1, p2, layout=(2, 1), size=(900, 600))
    end
    gif(anime, fps=8)
end
