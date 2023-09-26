#= 给定边界下，考虑存在恒定外力的多个小球运动：
1、小球完全相同，半径为R；
2、碰撞检测采用均匀网格法，网格边长为a，横向共Kx个网格，纵向共Ky个网格；
3、需要较小dt，避免粒子越界；
4、将测量体系的总机械能E，总动量Px、Py；
=#
using Plots, StatsBase, Combinatorics

# 主函数
function main()
    # 环境参数
    N = 10        # 总粒子数
    m = 1         # 粒子质量
    R = 5         # 粒子半径
    a = 2 * R     # 网格边长
    Kx = 30       # 横向网格数
    Ky = 10       # 纵向网格数
    Xmin = 0      # 边界
    Ymin = 0
    Xmax = Kx * a
    Ymax = Ky * a
    dt = 0.001    # 时间步长
    T = 50        # 最终时间

    # 初始条件
    x = Xmax * rand(N)
    y = Ymax * rand(N)
    u = 100 * rand(N) .- 50
    v = 100 * rand(N) .- 50
    f = 0.0            # 外力x向分量
    g = 0.0            # 外力y向分量

    # 显式步进============================================================================
    # 初始化：
    nmax = length(0:dt:T)    # 总时间层数
    uu = zeros(N)            # 预测的粒子速度
    vv = zeros(N)
    ufix = zeros(N)          # 经修正的粒子速度
    vfix = zeros(N)
    E = zeros(nmax)          # 各时间层总机械能
    U = zeros(N, nmax)       # 各时间层速度，位置信息
    V = zeros(N, nmax)
    Xp = zeros(N, nmax)
    Yp = zeros(N, nmax)

    U[:, 1] = u
    V[:, 1] = v
    Xp[:, 1] = x
    Yp[:, 1] = y
    E[1] = EPcalculate!(u, v, m, f, g, Xp, Yp, 1)

    # 更新 n >= 2  层的速度与位置：
    for n in 2:nmax
        uu = dt / m * f .+ U[:, n-1]   # 一阶时间精度
        vv = dt / m * g .+ V[:, n-1]
        x = dt * uu + Xp[:, n-1]
        y = dt * uu + Yp[:, n-1]
        ufix, vfix = CollisionDetection!(uu, vv, x, y, a, Kx, Xmin, Xmax, Ymin, Ymax, N, R)
        U[:, n] = ufix
        V[:, n] = vfix
        Xp[:, n] = dt * ufix + Xp[:, n-1]
        Yp[:, n] = dt * vfix + Yp[:, n-1]
        # 计算粒子群的机械能E，动量Px，Py
        E[n] = EPcalculate!(ufix, vfix, m, f, g, Xp, Yp, n)
    end

    return E, Xp, Yp, Xmin, Xmax, Ymin, Ymax, R, nmax
end

# 计算机械能，动量x分量，动量y分量
function EPcalculate!(ufix, vfix, m, f, g, Xp, Yp, n)
    E = sum(m / 2 * (ufix .^ 2 .+ vfix .^ 2)) + sum(g * (Yp[:, 1] - Yp[:, n]))
    return E
end

# 粒子间碰撞检测与修正
function CollisionDetection!(uu, vv, x, y, a, Kx, Xmin, Xmax, Ymin, Ymax, N, R)
    α = zeros(Int64, N)    # 所属网格的行序数
    β = zeros(Int64, N)    # 所属网格的列序数
    k = zeros(Int64, N)    # 所属网格的序数
    u = uu                 # 若无碰撞，则预测速度无需修正
    v = vv

    # 分析是否存在粒子碰壁，并修正：
    for i in 1:N
        if x[i] - Xmin <= R
            u[i] = abs(u[i])   # 采用abs()而非 -u[i] 为了让粒子远离边界，避免因速度过小而在下一时间步仍碰壁
        end
        if Xmax - x[i] <= R
            u[i] = -abs(u[i])
        end
        if y[i] - Ymin <= R
            v[i] = abs(v[i])
        end
        if Ymax - y[i] <= R
            v[i] = -abs(v[i])
        end
    end
    #

    # 遍历N个粒子，分析各粒子处于哪个网格内：
    for i in 1:N
        α[i] = Int64(ceil(x[i] / a))
        β[i] = Int64(ceil(y[i] / a))
        k[i] = Int64((α[i] - 1) * Kx + β[i])
    end

    # 挑选出包含2个及以上粒子的网格，分析是否存在粒子间碰撞：
    A = Array{Int64}[]
    B = Array{Int64}[]
    cp = StatsBase.countmap(k)      # 生成一个字典，统计k中每个值出现次数（键i：粒子所属格子的序数值k -> 值：粒子数）
    for i in keys(cp)
        if cp[i] >= 2
            A = findall(in(i), k)   # 找出满足条件的键在k中的索引(粒子索引)
            B = collect(combinations(A, 2))   # 从A中不重复地抽取2个索引的所有组合
            lb = length(B)
            for j in 1:lb
                r = sqrt((x[B[j][1]] - x[B[j][2]])^2 + (y[B[j][1]] - y[B[j][2]])^2)     # 球心距离
                if r <= 2R
                    u1, v1, u2, v2 = ParticleCollision01!(uu[B[j][1]], vv[B[j][1]], uu[B[j][2]], vv[B[j][2]], x[B[j][1]], y[B[j][1]], x[B[j][2]], y[B[j][2]])
                    u[B[j][1]] = u1
                    v[B[j][1]] = v1
                    u[B[j][2]] = u2
                    v[B[j][2]] = v2
                end
            end
        end
    end
    #

    return u, v
end

# 粒子间碰撞修正(考虑斜碰)：
function ParticleCollision01!(uu1, vv1, uu2, vv2, x1, y1, x2, y2)
    θ = atan((y2 - y1) / (x2 - x1))       # 球心连线与水平线夹角
    H = [cos(θ) sin(θ); -sin(θ) cos(θ)]   # 直角系向球心系幺正变换矩阵
    HT = H'                               # 逆变换矩阵
    ut1 = H * [uu1; vv1]
    ut2 = H * [uu2; vv2]
    u0 = copy(ut1)
    # 碰撞交换对心速度
    ut1[1] = ut2[1]
    ut2[1] = u0[1]
    # 速度逆变换
    uv1 = HT * ut1
    uv2 = HT * ut2
    u1 = uv1[1]
    v1 = uv1[2]
    u2 = uv2[1]
    v2 = uv2[2]

    return u1, v1, u2, v2
end

E, Xp, Yp, Xmin, Xmax, Ymin, Ymax, R, nmax = main();
# 绘制动画
begin
    Δn = 100
    p1 = Plots.plot(1, label="E/E₀", xlabel="steps / 100")
    δE = E ./ E[1]
    @gif for n in 1:Δn:nmax
        p2 = Plots.scatter(Xp[:, n], Yp[:, n],
            markersize=2R, markercolor=:red,
            aspect_ratio=1, legend=:false,
            xlim=[Xmin, Xmax], ylim=[Ymin, Ymax])

        push!(p1, δE[n])
        p3 = Plots.plot(p1, p2, layout=(2, 1), size=(800, 600))
    end
end
