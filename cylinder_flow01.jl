#= 采用LBGK模型，演算圆柱绕流问题
速度离散模型：D2Q9:
7---3---6
| \ | / |
4---1---2
| / | \ |
8---5---9
=#
using Plots

# 环境参数
H = 100                                    # 流域高SI[m/s]
W = 500                                    # 流域宽SI[m/s]
obstX = W / 5                              # 圆柱位置与半径
obstY = H / 2
obstR = H / 10
ρ₀ = 1.0                                   # 初始密度SI[kg/m³]
ν = 0.1                                    # 运动粘度SI[m²/s]
Re = 150                                   # 雷诺数
δx = 1                                     # 空间步长SI[m]
δt = 0.1                                   # 时间步长SI[s]
tmax = 1500                                # 计算时间SI[s]
u̅ᵢₙ = Re*ν/(2*obstR)                       # 入口平均速度SI[m/s]

# 求解域
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

# 存储设置
gap = 50                                                         # 每间隔gap个时间步存储一次数据用于返回
NUM = Int(round(M / gap)) + 1
save_ux = zeros(Ny + 1, Nx + 1, NUM)
save_uy = zeros(Ny + 1, Nx + 1, NUM)

# LBGK
function LBGK()
    # D2Q9:
    c₀ = δx / δt
    cₛ = c₀ / sqrt(3)
    τ = ν / (cₛ^2 * δt) + 1 / 2
    ω = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
    cx = c₀ * [0 1 0 -1 0 1 -1 -1 1]
    cy = c₀ * [0 0 1 0 -1 1 1 -1 -1]

    # 初始化:
    ρ = ones(Ny + 1, Nx + 1) * ρ₀
    u = zeros(Ny + 1, Nx + 1, 2)
    U = 6 * u̅ᵢₙ / H^2 * Y[2:Ny, 1] .* (H .- Y[2:Ny, 1])           # 入口速度抛物型分布设置
    u[2:Ny, 1, 1] = U
    f = zeros(Ny + 1, Nx + 1, 9)
    F = zeros(Ny + 1, Nx + 1, 9)
    # 分布函数初值:
    for k in 1:9
        f[:, :, k] = ω[k] * ρ[:, :] .* (1 .+ 1 / cₛ^2 * (
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) +
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) .^ 2 * (1 / (2 * cₛ^2)) -
            (u[:, :, 1] .^ 2 + u[:, :, 2] .^ 2) * (1 / 2)
        ))
    end

    # 定义平衡态函数:
    function fEq1(rho::Matrix, u::Array{Float64, 3}, k::Int)
        feq = ω[k] * rho .* (1 .+ 1 / cₛ^2 * (
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) +
            (cx[k] * u[:, :, 1] + cy[k] * u[:, :, 2]) .^ 2 * (1 / (2 * cₛ^2)) -
            (u[:, :, 1] .^ 2 + u[:, :, 2] .^ 2) * (1 / 2)
        ))
        return feq
    end
    function fEq2(rho::Vector, u::Vector, v::Vector, k::Int)
        feq = ω[k] * rho .* (1 .+ 1 / cₛ^2 * (
            (cx[k] * u + cy[k] * v) +
            (cx[k] * u + cy[k] * v) .^ 2 * (1 / (2 * cₛ^2)) -
            (u .^ 2 + v .^ 2) * (1 / 2)
        ))
        return feq
    end

    # LBGK:
    num = 1    # 存储计数
    for tk in 1:M
        # 碰撞与迁移(内点)
        f[2:Ny, 2:Nx, 1] = f[2:Ny, 2:Nx, 1] + 1 / τ * (fEq1(ρ[2:Ny, 2:Nx], u[2:Ny, 2:Nx, :], 1) - f[2:Ny, 2:Nx, 1])
        f[2:Ny, 2:Nx, 2] = f[2:Ny, 1:Nx-1, 2] + 1 / τ * (fEq1(ρ[2:Ny, 1:Nx-1], u[2:Ny, 1:Nx-1, :], 2) - f[2:Ny, 1:Nx-1, 2])
        f[2:Ny, 2:Nx, 3] = f[1:Ny-1, 2:Nx, 3] + 1 / τ * (fEq1(ρ[1:Ny-1, 2:Nx], u[1:Ny-1, 2:Nx, :], 3) - f[1:Ny-1, 2:Nx, 3])
        f[2:Ny, 2:Nx, 4] = f[2:Ny, 3:Nx+1, 4] + 1 / τ * (fEq1(ρ[2:Ny, 3:Nx+1], u[2:Ny, 3:Nx+1, :], 4) - f[2:Ny, 3:Nx+1, 4])
        f[2:Ny, 2:Nx, 5] = f[3:Ny+1, 2:Nx, 5] + 1 / τ * (fEq1(ρ[3:Ny+1, 2:Nx], u[3:Ny+1, 2:Nx, :], 5) - f[3:Ny+1, 2:Nx, 5])
        f[2:Ny, 2:Nx, 6] = f[1:Ny-1, 1:Nx-1, 6] + 1 / τ * (fEq1(ρ[1:Ny-1, 1:Nx-1], u[1:Ny-1, 1:Nx-1, :], 6) - f[1:Ny-1, 1:Nx-1, 6])
        f[2:Ny, 2:Nx, 7] = f[1:Ny-1, 3:Nx+1, 7] + 1 / τ * (fEq1(ρ[1:Ny-1, 3:Nx+1], u[1:Ny-1, 3:Nx+1, :], 7) - f[1:Ny-1, 3:Nx+1, 7])
        f[2:Ny, 2:Nx, 8] = f[3:Ny+1, 3:Nx+1, 8] + 1 / τ * (fEq1(ρ[3:Ny+1, 3:Nx+1], u[3:Ny+1, 3:Nx+1, :], 8) - f[3:Ny+1, 3:Nx+1, 8])
        f[2:Ny, 2:Nx, 9] = f[3:Ny+1, 1:Nx-1, 9] + 1 / τ * (fEq1(ρ[3:Ny+1, 1:Nx-1], u[3:Ny+1, 1:Nx-1, :], 9) - f[3:Ny+1, 1:Nx-1, 9])
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
        # 宏观量ρ 、u
        ρ = sum(f, dims=3)
        u[:, :, 1] = sum(reshape(cx .* reshape(f, (Ny + 1) * (Nx + 1), 9), Ny + 1, Nx + 1, 9), dims=3) ./ ρ
        u[:, :, 2] = sum(reshape(cy .* reshape(f, (Ny + 1) * (Nx + 1), 9), Ny + 1, Nx + 1, 9), dims=3) ./ ρ
        # 入口，出口的非平衡外推边界(要基于边界与近邻流体格点的新分布值)
        # 出口x向压力梯度为0，故出口速度x向梯度为0
        ρ[:, Nx+1] = ρ[:, Nx]
        u[2:Ny, 1, 1] = U
        u[2:Ny, 1, 2] .= 0.0
        u[:, Nx+1, 1] = u[:, Nx, 1]
        u[:, Nx+1, 2] .= 0.0
        # 边界分布的非平衡外推(Guo格式)
        for k in 1:9
            f[2:Ny, 1, k] = fEq2(ρ[2:Ny, 1], u[2:Ny, 1, 1], u[2:Ny, 1, 2], k) + f[2:Ny, 2, k] - fEq2(ρ[2:Ny, 2], u[2:Ny, 2, 1], u[2:Ny, 2, 2], k)
            f[2:Ny, Nx+1, k] = fEq2(ρ[2:Ny, Nx+1], u[2:Ny, Nx+1, 1], u[2:Ny, Nx+1, 2], k) + f[2:Ny, Nx, k] - fEq2(ρ[2:Ny, Nx], u[2:Ny, Nx, 1], u[2:Ny, Nx, 2], k)
        end
        # 保存各向速度ux，uy
        if tk == (num - 1) * gap + 1
            save_ux[:, :, num] = u[:, :, 1]
            save_uy[:, :, num] = u[:, :, 2]
            num += 1
        end
    end

    return save_ux, save_uy
end
save_ux, save_uy = @time LBGK()

# 绘图:
umod = zeros(Ny + 1, Nx + 1, NUM)
Ω = zeros(Ny+1, Nx+1, NUM)
# 速度模SI[m/s]
for num in 1:NUM
    umod[:, :, num] = sqrt.(save_ux[:, :, num] .^ 2 + save_uy[:, :, num] .^ 2)
end
# 涡量场SI[1/s]
for num in 1:NUM
    for i in 2:Ny
        for j in 2:Nx
            Ω[i, j, num] = 1 / (2 * δx) * ((save_uy[i, j+1, num] - save_uy[i, j-1, num]) - (save_ux[i+1, j, num] - save_ux[i-1, j, num]))
        end
    end
    Ω[:, 1, num] = Ω[:, 2, num]
    Ω[:, Nx+1, num] = Ω[:, Nx, num]
    Ω[1, :, num] = Ω[2, :, num]
    Ω[Ny+1, :, num] = Ω[Ny, :, num]
end
# 动态图(Plots)
begin
    anime = @animate for num in 1:NUM
        ttl = string("Re = ", Re, "\n","t = ", (num-1)*gap*δt, "s", "\n", "Speed")
        p1 = Plots.heatmap(x, y, umod[:, :, num],
            color=:turbo,
            title=ttl,
            aspect_ratio=1,
            xlims=(0, W), ylims=(0, H), 
            colorbar_title = "m/s",
            xlabel="x (m)", ylabel="y (m)")
            
        
        p2 = Plots.heatmap(x, y, Ω[:, :, num],
            color=:berlin,
            aspect_ratio=1,
            title="Vorticity field",
            xlims=(0, W), ylims=(0, H), 
            clims=(-0.2, 0.2),
            colorbar_title = "1/s",
            xlabel="x (m)", ylabel="y (m)")
        
        p3 = Plots.plot(p1, p2, layout=(2, 1), size=(900, 600))
    end
    gif(anime, fps=8)
end
