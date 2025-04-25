using GLMakie
using Observables
using Printf
using LinearAlgebra
using GeometryBasics: Point2f


#=function wave_speed7(x, y)
    r = sqrt(x^2 + y^2) / 2
    if 7.45 <= r <= 10.0
        return (11.0 - r)/10
    elseif 6.5 <= r < 7.45
        return (160.0 - 21.0 * r)/10
    elseif 0.0 <= r < 6.5
        return (30.0 - r)/10
    else
        return 0.0
    end
end=#

#=function wave_speed2(x,y)
    r = norm([x,y])
    if r <= 10
        return 3.0
    else
        return 0.0
    end
end=#

# define wave velocity distributions
function wave_speed1(x,y) # homogenous
    if y > 0
        return 0
    else
        return 3
    end
end

function wave_speed2(x,y) # 2 layer
    if y > 0
        return 0
    elseif y > -10
        return 1.0
    else
        return 3.0
    end
end

function wave_speed3(x,y) # linear
    if y > 0
        return 0
    else
        return 2-y
    end
end

function wave_speed4(x,y) # low verocity zone
    if y > 0
        return 0
    else
        return 2-y - 5.0*exp(-((x-15)^2)/75 - ((y+10)^2)/25)
    end
end

function wave_speed5(x,y) # velocity up, down, and up as go deep
    if y > 0
        return 0
    elseif y > -4
        return (1 - y)*2.0
    elseif y > -8
        return (9 + y)*2.0
    else
        return (-y - 7)*2.0
    end
end

function wave_speed6(x,y) # abnormal seismic intensity(異常震域)
    if y > 0
        return 0
    else
        return 2 + 4*exp((-(x-y-40)^2)/100)
    end
end

# define gradients
function gradv_x(x, y, speed_func)
    dx = 0.00001
    return (speed_func(x+dx,y)-speed_func(x-dx,y))/(2*dx)
end

function gradv_y(x, y, speed_func)
    dy = 0.00001
    return (speed_func(x,y+dy)-speed_func(x,y-dy))/(2*dy)
end


# create velocity distribution map
function speedmap(n_points, speed_func)
    dx = (0.0 - (-100.0)) / n_points
    dy = (20.0 - (-5.0)) / n_points

    speed_x = zeros(n_points)
    speed_y = zeros(n_points)
    speed_val = zeros(n_points,n_points)

    for i in 1:n_points
        speed_x[i] =  dx * i
        speed_y[i] = -20.0 + dy * i
    end

    for i in 1:n_points
        for j in 1:n_points
            speed_val[i, j] = speed_func(speed_x[i], speed_y[j])
        end
    end
    return speed_val
end

# change elements<0.01 to NaN
function process_matrix(data)
    # 値が 0 の部分を NaN にすることで、マッピングしない
    masked = map(x -> x < 0.01 ? NaN : x, data)
    return masked
end

# define ray equation
function rayequation(x, y, theta, dt, speed_func_obs)
    speed_func = speed_func_obs[]

    x_new = x .+ speed_func.(x, y) .* cos.(theta) .* dt
    y_new = y .+ speed_func.(x, y) .* sin.(theta) .* dt

    theta_new = theta .+ dt .* (
        gradv_x.(x_new, y_new, speed_func) .* sin.(theta) .-
        gradv_y.(x_new, y_new, speed_func) .* cos.(theta)
    )

    return x_new, y_new, theta_new
end

# main function to set window
function showwindow()

    # define Constants and Variables
    n_points = 1000 
    n_angles = 10001 #wave points +1

    t = Observable(0.00)
    dt = 0.01

    x = Observable(fill(2.0, n_angles))
    y = Observable(fill(-2.0, n_angles))
    x_int = Observable(fill(2.0, n_angles))
    y_int = Observable(fill(-2.0, n_angles))

    x0 = Observable(2.0)
    y0 = Observable(2.0)

    theta = zeros(n_angles)
    theta_int = zeros(n_angles)

    t_reach = Observable(fill(NaN, n_angles))
    x_reach = Observable(fill(NaN, n_angles))

    # incident angles
    for i in 1:n_angles
        theta[i] = 2π * (i - 1) / (n_angles-1)
        theta_int[i] = 2π * (i - 1) / (n_angles-1)
    end

    fig = Figure(
        size=(1000,600),
        backgroundcolor = :white
    )

    ax1 = Axis(
        fig[1:9,1:10],
        aspect = AxisAspect(4),
        xrectzoom=false, yrectzoom=false,
        limits = (0,88,-20,2)
    )

    ax2 = Axis(
        fig[11:19,1:10],
        aspect = AxisAspect(4),
        xrectzoom=false, yrectzoom=false,
        limits = (0,88,0,11)
    )

    # travel time plot
    ttpoints = Observable(Point2f.(x_reach[], t_reach[]))
    ttplot = scatter!(ax2,ttpoints, markersize=1.0, color="red")
    #ttplot_theta_data = Observable(t_reach[]) # ← to confirm the accuracy
    #ttplot_theta = scatter!(ax2, 1:n_angles, ttplot_theta_data, color=:red, markersize=1.0) # ← to confirm the accuracy


    # Start/Stop button
    start = Observable(false)
    bt_startstop = Button(fig[2,11], label="start")
    on(bt_startstop.clicks) do _
        if start[]
            start[] = false
            bt_startstop.label = "start"
            #println(t_reach[]) # ← to confirm the accuracy
        else
            start[] = true
            bt_startstop.label = "stop"
        end
    end

    # time text
    timerlabel = Label(
        fig[3,11],
        text = lift(x -> @sprintf("%.2fs", x), t),
        fontsize = 30,
        halign = :center
    )
    
    # reset button
    bt_reset = Button(fig[4,11],label="reset")
    on(bt_reset.clicks) do _
        x[] .= x_int[]
        y[] .= y_int[]
        theta .= theta_int
        t[] = 0.0
        raypoints[] = Point2f.(x[], y[])
        x_reach[] .= fill(NaN, n_angles)
        t_reach[] .= fill(NaN, n_angles)
        ttpoints[] = Point2f.(x_reach[], t_reach[])
    end

    # set initial position (hypocenter)
    Label(
        fig[5,11],
        text = "震源のx座標",
        fontsize=20,
        font = "Hiragino Kaku Gothic ProN"
    )
    x0_slider = Slider(fig[6,11], range = 0.0:1.0:80.0, startvalue = 2.0)
    on(x0_slider.value) do val
        x0[] = val
        x[] .= fill(x0[], n_angles)
        x_int[] .= fill(x0[], n_angles)
        initial_points[] = Point2f.(x_int[], y_int[])
    end

    Label(
        fig[7,11],
        text = "震源のy座標",
        fontsize=20,
        font = "Hiragino Kaku Gothic ProN"
    )
    y0_slider = Slider(fig[8,11], range = 0.0:0.5:20.0, startvalue = 2.0)
    on(y0_slider.value) do val
        y0[] = -val
        y[] .= fill(y0[], n_angles)
        y_int[] .= fill(y0[], n_angles)
        initial_points[] = Point2f.(x_int[], y_int[])
    end


    # title text
    Label(fig[0,1], "地震波シミュレーション", fontsize=40, halign=:left, font = "Hiragino Kaku Gothic ProN")

    # wave velocity selection
    wave_speed_func = Observable{Function}(wave_speed1)
    initial_matrix = process_matrix(speedmap(n_points, wave_speed_func[]))
    hm = Observable(heatmap!(ax1, LinRange(0.0, 100.0, n_points), LinRange(-20.0, 5.0, n_points), initial_matrix, colormap = :cool))
    Colorbar(fig[9, 2:7], hm[], label = "Wave Speed", vertical = false)

    selected_label = Label(fig[2,12], "現在の波速:速度分布1", fontsize=20, font = "Hiragino Kaku Gothic ProN")

    wave_speed_funcs = [wave_speed1, wave_speed2, wave_speed3, wave_speed4, wave_speed5, wave_speed6]
    labels = ["速度一様", "水平2層", "水平成層", "速度分布4", "速度分布5","速度分布6"]

    for i in 1:6
        btn = Button(fig[2+i, 12], label=labels[i], font = "Hiragino Kaku Gothic ProN")
        on(btn.clicks) do _
            wave_speed_func[] = wave_speed_funcs[i]
            selected_label.text[] = "現在の波速:" * labels[i]
            new_map = process_matrix(speedmap(n_points, wave_speed_func[]))
            hm[].values = new_map
        end
    end

    # plot hypocenter
    initial_points = Observable(Point2f.(x_int[], y_int[]))
    initial_position = scatter!(ax1, initial_points, markersize=10, color="black")

        
    # show hypocenter x position to ax2
    vline_pos = Observable(x0[])
    vline = vlines!(ax2, vline_pos, color=:blue, linewidth=2.0)
    on(x0_slider.value) do val
        vline_pos[] = val
    end

    # plot wave points
    raypoints = Observable(Point2f.(x[], y[]))
    scatterplot = scatter!(ax1,raypoints,markersize=1.0, color="black")

    display(fig)

    # apply ray equation
    while true
        if start[]
            sleep(dt)
            t[] += dt
            x[],y[],theta = rayequation(x[],y[],theta,dt,wave_speed_func)
            raypoints[] = Point2f.(x[], y[])
            for i in 1:n_angles
                if isnan(x_reach[][i]) && y[][i] > 0
                    x_reach[][i] = x[][i]
                    t_reach[][i] = t[]
                end
            end
            ttpoints[] = Point2f.(x_reach[], t_reach[])
            #ttplot_theta_data[] = copy(t_reach[]) # ← to confirm the accuracy
        else
            sleep(0.1)
        end
    end

end

showwindow()
