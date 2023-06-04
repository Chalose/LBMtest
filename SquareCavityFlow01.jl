#= 参考何雅玲《格子Boltzmann方法的理论及应用》附录D顶盖驱动代码(LBGK)
D2Q9模型:
7---3---6
| \ | / |
4---1---2
| / | \ |
8---5---9
边界采用非平衡外推(Guo)
绘图采用Syslab提供的TyPlot包
=#

# 环境参数(各参数均为格子单位，即δx = 1, δt = 1)===================================================
Nx = 100                    # x向格子数
Ny = 100                    # y向格子数
U = 0.1                     # 顶盖速度
Re = 1000                   # 雷诺数
ν = U * Nx / Re             # 动力粘度
ρ₀ = 1.0                    # 初始密度
tmax = 40000                # 最大求解时间

# 求解域==========================================================================================
x = LinRange(0, Nx + 1, Nx + 1)
y = LinRange(0, Ny + 1, Ny + 1)
X = x' .* ones(Ny + 1, Nx + 1)           # 右为x正向，下为y正向
Y = ones(Ny + 1, Nx + 1) .* y
t = 0:1:tmax
M = length(t)

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

# 主函数===========================================================================================
function mainFun()
    # 初始化
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
    @views for k in 1:M
        # 内点碰撞与迁移............................................................................
        f[2:Ny, 2:Nx, 1] = f[2:Ny, 2:Nx, 1] * (1 - 1 / τ) + 1 / τ * fEq!(u[2:Ny, 2:Nx], v[2:Ny, 2:Nx], ρ[2:Ny, 2:Nx], 1)
        f[2:Ny, 2:Nx, 2] = f[2:Ny, 1:Nx-1, 2] * (1 - 1 / τ) + 1 / τ * fEq!(u[2:Ny, 1:Nx-1], v[2:Ny, 1:Nx-1], ρ[2:Ny, 1:Nx-1], 2)
        f[2:Ny, 2:Nx, 3] = f[1:Ny-1, 2:Nx, 3] * (1 - 1 / τ) + 1 / τ * fEq!(u[1:Ny-1, 2:Nx], v[1:Ny-1, 2:Nx], ρ[1:Ny-1, 2:Nx], 3)
        f[2:Ny, 2:Nx, 4] = f[2:Ny, 3:Nx+1, 4] * (1 - 1 / τ) + 1 / τ * fEq!(u[2:Ny, 3:Nx+1], v[2:Ny, 3:Nx+1], ρ[2:Ny, 3:Nx+1], 4)
        f[2:Ny, 2:Nx, 5] = f[3:Ny+1, 2:Nx, 5] * (1 - 1 / τ) + 1 / τ * fEq!(u[3:Ny+1, 2:Nx], v[3:Ny+1, 2:Nx], ρ[3:Ny+1, 2:Nx], 5)
        f[2:Ny, 2:Nx, 6] = f[1:Ny-1, 1:Nx-1, 6] * (1 - 1 / τ) + 1 / τ * fEq!(u[1:Ny-1, 1:Nx-1], v[1:Ny-1, 1:Nx-1], ρ[1:Ny-1, 1:Nx-1], 6)
        f[2:Ny, 2:Nx, 7] = f[1:Ny-1, 3:Nx+1, 7] * (1 - 1 / τ) + 1 / τ * fEq!(u[1:Ny-1, 3:Nx+1], v[1:Ny-1, 3:Nx+1], ρ[1:Ny-1, 3:Nx+1], 7)
        f[2:Ny, 2:Nx, 8] = f[3:Ny+1, 3:Nx+1, 8] * (1 - 1 / τ) + 1 / τ * fEq!(u[3:Ny+1, 3:Nx+1], v[3:Ny+1, 3:Nx+1], ρ[3:Ny+1, 3:Nx+1], 8)
        f[2:Ny, 2:Nx, 9] = f[3:Ny+1, 1:Nx-1, 9] * (1 - 1 / τ) + 1 / τ * fEq!(u[3:Ny+1, 1:Nx-1], v[3:Ny+1, 1:Nx-1], ρ[3:Ny+1, 1:Nx-1], 9)
    
        # 内点宏观量ρ, u, v刷新......................................................................
        @views begin
            ρ[2:Ny, 2:Nx] = sum!(ρ[2:Ny, 2:Nx], f[2:Ny, 2:Nx, :])
            u[2:Ny, 2:Nx] = sum!(u[2:Ny, 2:Nx], reshape(cx .* reshape(f[2:Ny, 2:Nx, :], (Ny - 1) * (Nx - 1), 9), Ny - 1, Nx - 1, 9)) ./ ρ[2:Ny, 2:Nx]
            v[2:Ny, 2:Nx] = sum!(v[2:Ny, 2:Nx], reshape(cy .* reshape(f[2:Ny, 2:Nx, :], (Ny - 1) * (Nx - 1), 9), Ny - 1, Nx - 1, 9)) ./ ρ[2:Ny, 2:Nx]
        end
    
        # 边界处理(非平衡外推).......................................................................
        # 确定边界宏观量
        @views begin
            ρ[:, 1] = ρ[:, 2]          # left
            ρ[:, Nx+1] = ρ[:, Nx]      # right
            ρ[1, :] = ρ[2, :]          # bottom
            ρ[Ny+1, :] = ρ[Ny, :]      # top
    
            #u[:, [1, Nx + 1]] .= 0.0  # 初始u, v零矩阵，可省略
            #v[:, [1, Nx + 1]] .= 0.0
            #u[1, :] .= 0.0
            #v[1, :] .= 0.0
            u[Ny+1, :] .= U
            #v[Ny+1, :] .= 0.0
        end
        # 边界分布非平衡外推(Guo格式)
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
    end

    return u, v
end
u, v = @time mainFun();

# 绘图=============================================================================================
uv = @. sqrt(u^2 + v^2)
# TyPlot
begin
    subplot(121)
    p1 = TyPlot.pcolor(x, y, uv)
    TyPlot.axis("equal", "tight")
    TyPlot.colormap(p1, "jet")
    #TyPlot.colorbar(p1)
    tit1 = string("speed", '\n', "Re = ", Re)
    title(tit1)

    subplot(122)
    startx1 = 80:0.01:Nx+1    # 右下角涡
    starty1 = 0:0.01:21
    TyPlot.streamline(X, Y, u, v, startx1, starty1, color="black")
    hold("on")
    startx2 = 1:0.05:Nx+1     # 大涡与左下角涡
    starty2 = 1:0.05:Ny+1
    TyPlot.streamline(X, Y, u, v, startx2, starty2, color="black")
    Xx, Yy = meshgrid2(2:5:Nx+1, 2:5:Ny+1)
    TyPlot.quiver(Xx, Yy, u[2:5:Nx+1, 2:5:Ny+1], v[2:5:Nx+1, 2:5:Ny+1], length=200)
    hold("off")
    TyPlot.axis("equal", "tight")
    tit2 = string("streamline", '\n', "Re = ", Re)
    title(tit2)
end
