using Random
using LinearAlgebra
using SpecialFunctions
using Dates
using Base: Float64, Integer
using Distributions
using Printf
using CairoMakie
using LaTeXStrings
using Statistics
au_to_ev = 27.211


##================================================================##
function marginals_helper(mat)
##================================================================##
    nx, ny  = size(mat)
    marg_x  = zeros(nx)
    marg_ir = zeros(ny)
    for i=1:nx
        marg_x[i]  = sum(mat[i, :])
        marg_ir[i] = sum(mat[:,i])
    end

    return marg_x, marg_ir
end
##================================================================##
function plt_signal_decomposition_v2(omega_x::Vector{Float64}, omega_ir::Vector{Float64}, 
                            mat_incoh::Matrix{Float64}, mat_coh::Matrix{Float64}, outfile_1::String, outfile_2::String)
##================================================================##
    marg_x_incoh, marg_ir_incoh = marginals_helper(mat_incoh)
    idx_x  = sortperm(marg_x_incoh, rev=true)
    idx_ir = sortperm(marg_ir_incoh, rev=true)
    omega_x  = omega_x .* au_to_ev
    omega_ir = omega_ir .* au_to_ev
    n::Int = 20
    for islice=1:n
        ix = idx_x[islice]
        ir = idx_ir[islice]

        if marg_ir_incoh[ir] <= 1.0e-3
            continue
        else
            val_x    = round(Int, omega_x[ix])
            val_ir   = round(omega_ir[ir]; digits=4)


            fig_1 = Figure(size = (1200, 900))
            ax1 = Axis(fig_1[1,1], 
                xlabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}", ylabel = LaTeXString("P(\\omega_{\\mathrm{X}}=$(val_x)\\mathrm{eV},\\,\\omega_{\\mathrm{IR}})"),
                xlabelsize = 36, ylabelsize = 36,
                xgridvisible = true, ygridvisible = true,
                xminorgridvisible = false, yminorgridvisible = false,
                xticksvisible = true, yticksvisible = true,
                xticklabelsize = 24, yticklabelsize = 24,
                xminorticksvisible = true, yminorticksvisible = true,
                xticksmirrored = true, yticksmirrored = true,
                xtickalign = 0.66, ytickalign = 0.66,
                xminortickalign = 1, yminortickalign = 1,
                xticksize = 15, yticksize = 15,
                xminorticksize = 5, yminorticksize = 5, 
                xlabelpadding = 18, ylabelpadding = 18,
                )
            

            interfer_vec = mat_coh[ix, :] - mat_incoh[ix, :]
            lines!(ax1, omega_ir, mat_incoh[ix, :], linewidth=2.0, color = :darkred, label = "Incoherent")
            lines!(ax1, omega_ir, mat_coh[ix, :], linewidth=2.0, color = :darkblue, label = "Coherent")
            lines!(ax1, omega_ir, interfer_vec, linewidth=2.0, color = :orange, label = "Interference")

            leg = Legend(fig_1[1,1], ax1, patchsize = (30,20), padding = (15, 15, 15, 15), margin = (15,15,15,15), 
                halign = :right, valign = :top, tellwidth = false,
                labelsize = 28, markersize = 26)
               
            outfile_ir = outfile_1 * "_$(islice).pdf"
            save(outfile_ir, fig_1, px_per_unit = 2)
        end

        #-------
        if marg_x_incoh[ix] <= 1.0e-3
            continue
        else
            fig_2 = Figure(size = (1200, 900))
            ax2 = Axis(fig_2[1,1], 
                xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", ylabel = LaTeXString("P(\\omega_{\\mathrm{X}}\\mathrm{eV},\\,\\omega_{\\mathrm{IR}}=$(val_ir))"),
                xlabelsize = 36, ylabelsize = 36,
                xgridvisible = true, ygridvisible = true,
                xminorgridvisible = false, yminorgridvisible = false,
                xticksvisible = true, yticksvisible = true,
                xticklabelsize = 24, yticklabelsize = 24,
                xminorticksvisible = true, yminorticksvisible = true,
                xticksmirrored = true, yticksmirrored = true,
                xtickalign = 0.66, ytickalign = 0.66,
                xminortickalign = 1, yminortickalign = 1,
                xticksize = 15, yticksize = 15,
                xminorticksize = 5, yminorticksize = 5, 
                xlabelpadding = 18, ylabelpadding = 18,
                )

            # xlim = 45:75
            # omega_x = omega_x[xlim]
            # matrix  = matrix[:, ir][xlim]
            # matrix = matrix[:, ir+1]

            interfer_vec = mat_coh[:, ir] - mat_incoh[:, ir]
            lines!(ax2, omega_x, mat_incoh[:, ir], linewidth=3.0, color = :darkred, label = "Incoherent")
            lines!(ax2, omega_x, mat_coh[:, ir], linewidth=2.8, color = :darkblue, label = "Coherent")
            lines!(ax2, omega_x, interfer_vec, linewidth=3.0, color = :orange, label = "Interference")

            leg = Legend(fig_2[1,1], ax2, patchsize = (30,20), padding = (15, 15, 15, 15), margin = (15,15,15,15), 
                halign = :right, valign = :top, tellwidth = false,
                labelsize = 28, markersize = 26)
                    
            outfile_x = outfile_2 * "_$(islice).pdf"
            save(outfile_x, fig_2, px_per_unit = 2)
        end
    end
    return nothing
end
##================================================================##
function plt_signal_decomposition_v3(omega_x::Vector{Float64}, omega_ir::Vector{Float64}, 
                            mat_incoh_zo::Matrix{Float64}, mat_coh_zo::Matrix{Float64}, mat_incoh_fo::Matrix{Float64}, mat_coh_fo::Matrix{Float64}, outfile_1::String)
##================================================================##
    max_val, pair = findmax(mat_incoh_zo)
    ix = pair[1]
    ir = pair[2]

    omega_x  = omega_x .* au_to_ev
    omega_ir = omega_ir .* au_to_ev
    val_x    = round(Int, omega_x[ix])
    val_ir   = round(omega_ir[ir]; digits=4)


    fig_1 = Figure(size = (1200, 900))
    ax1   = Axis(fig_1[1,1], 
            xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", 
            ylabel = LaTeXString("P(\\omega_{\\mathrm{X}}\\,\\omega_{\\mathrm{IR}}=$(val_ir) \\mathrm{eV},)"),
            xlabelsize = 36, ylabelsize = 36,
            xgridvisible = true, ygridvisible = true,
            xminorgridvisible = false, yminorgridvisible = false,
            xticksvisible = true, yticksvisible = true,
            xticklabelsize = 24, yticklabelsize = 24,
            xminorticksvisible = true, yminorticksvisible = true,
            xticksmirrored = true, yticksmirrored = true,
            xtickalign = 0.66, ytickalign = 0.66,
            xminortickalign = 1, yminortickalign = 1,
            xticksize = 15, yticksize = 15,
            xminorticksize = 5, yminorticksize = 5, 
            xlabelpadding = 18, ylabelpadding = 18,
            )
    xlim = 45:75
    omega_x = omega_x[xlim]
    colors = (
    incoh_first = "#1F4E79",  # deep steel blue
    coh_first   = "#8C2D2D",  # muted dark red

    incoh_zero  = "#7FA6C9",  # soft desaturated blue
    coh_zero    = "#C98282"   # soft desaturated red
)
    lines!(ax1, omega_x, mat_incoh_zo[xlim, ir], linewidth=3.2, color = colors.incoh_zero, label = "Zero-order Incoherent")
    lines!(ax1, omega_x, mat_coh_zo[xlim, ir], linewidth=2.8, color = colors.coh_zero, label = "Zero-order Coherent")
    lines!(ax1, omega_x, mat_incoh_fo[xlim, ir], linewidth=3.2, color = colors.incoh_first, label = "First-order Incoherent")
    lines!(ax1, omega_x, mat_coh_fo[xlim, ir], linewidth=2.8, color = colors.coh_first, label = "First-order Coherent")

    leg = Legend(fig_1[1,1], ax1, patchsize = (30,20), padding = (15, 15, 15, 15), margin = (15,15,15,15), 
        halign = :right, valign = :top, tellwidth = false,
        labelsize = 28, markersize = 26)

    save(outfile_1, fig_1, px_per_unit = 2)
    return nothing
end
##================================================================##
function plt_signal_decomposition_v1(omega_x::Vector{Float64}, omega_ir::Vector{Float64}, 
                            mat_incoh::Matrix{Float64}, mat_coh::Matrix{Float64}, outfile_1::String, outfile_2::String)
##================================================================##
    max_val, pair = findmax(mat_incoh)
    ix = pair[1]
    ir = pair[2]

    omega_x  = omega_x .* au_to_ev
    omega_ir = omega_ir .* au_to_ev
    val_x    = round(Int, omega_x[ix])
    val_ir   = round(omega_ir[ir]; digits=4)


    fig_1 = Figure(size = (1200, 900))
    ax1   = Axis(fig_1[1,1], 
            xlabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}", 
            ylabel = LaTeXString("P(\\omega_{\\mathrm{X}}=$(val_x) \\mathrm{eV},\\,\\omega_{\\mathrm{IR}})"),
            xlabelsize = 36, ylabelsize = 36,
            xgridvisible = true, ygridvisible = true,
            xminorgridvisible = false, yminorgridvisible = false,
            xticksvisible = true, yticksvisible = true,
            xticklabelsize = 24, yticklabelsize = 24,
            xminorticksvisible = true, yminorticksvisible = true,
            xticksmirrored = true, yticksmirrored = true,
            xtickalign = 0.66, ytickalign = 0.66,
            xminortickalign = 1, yminortickalign = 1,
            xticksize = 15, yticksize = 15,
            xminorticksize = 5, yminorticksize = 5, 
            xlabelpadding = 18, ylabelpadding = 18,
            )
    

    interfer_vec = mat_coh[ix, :] - mat_incoh[ix, :]
    lines!(ax1, omega_ir, mat_incoh[ix, :], linewidth=3.2, color = :darkred, label = "Incoherent")
    lines!(ax1, omega_ir, mat_coh[ix, :], linewidth=2.8, color = :darkblue, label = "Coherent")
    lines!(ax1, omega_ir, interfer_vec, linewidth=3.2, linestyle=(:dot, :dense), color = :orange, label = "Interference")

    leg = Legend(fig_1[1,1], ax1, patchsize = (30,20), padding = (15, 15, 15, 15), margin = (15,15,15,15), 
        halign = :left, valign = :top, tellwidth = false,
        labelsize = 28, markersize = 26)

    save(outfile_1, fig_1, px_per_unit = 2)

    #-------
    # fig_2 = Figure(size = (1200, 900))
    # ax2   = Axis(fig_2[1,1], 
    #     xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", 
    #     ylabel = LaTeXString("P(\\omega_{\\mathrm{X}},\\,\\omega_{\\mathrm{IR}}=$(val_ir) \\mathrm{eV})"),
        # xlabelsize = 36, ylabelsize = 36,
        # xgridvisible = true, ygridvisible = true,
        # xminorgridvisible = false, yminorgridvisible = false,
        # xticksvisible = true, yticksvisible = true,
        # xticklabelsize = 24, yticklabelsize = 24,
        # xminorticksvisible = true, yminorticksvisible = true,
        # xticksmirrored = true, yticksmirrored = true,
        # xtickalign = 0.66, ytickalign = 0.66,
        # xminortickalign = 1, yminortickalign = 1,
        # xticksize = 15, yticksize = 15,
        # xminorticksize = 5, yminorticksize = 5, 
        # xlabelpadding = 18, ylabelpadding = 18,
        # )

    # # xlim = 45:75
    # # omega_x = omega_x[xlim]
    # # matrix  = matrix[:, ir][xlim]
    # # matrix = matrix[:, ir+1]

    # interfer_vec = mat_coh[:, ir] - mat_incoh[:, ir]
    # lines!(ax2, omega_x, mat_incoh[:, ir], linewidth=2.0, color = :darkred, label = "Incoherent")
    # lines!(ax2, omega_x, mat_coh[:, ir], linewidth=2.0, color = :darkblue, label = "Coherent")
    # lines!(ax2, omega_x, interfer_vec, linewidth=2.0, color = :orange, label = "Interference")

    # leg = Legend(fig_2[1,1], ax2, patchsize = (20,20), padding = (10,10,10,10), margin = (10,10,10,10), 
    #     halign = :left, valign = :top, tellwidth = false,
    #     labelsize = 20, markersize = 20)
    # save(outfile_2, fig_2)

    return nothing
end
##================================================================##
function get_ir_region(ir_val)
##================================================================##
    if 0.886 <= ir_val <= 1.700
        return "NIR"
    elseif 0.413 <= ir_val <= 0.8859
        return "SWIR"
    elseif 0.155 <= ir_val <= 0.4129
        return "MWIR"
    elseif 0.083 <= ir_val <= 0.1549
        return "LWIR"
    elseif 0.0012 <= ir_val <= 0.0829
        return "FIR"
    else
        println("$ir_val is not in IR range!")
        return nothing
    end
end
##================================================================##
function heatmap_plot_matrix_representation_v2_gamma(
    au_to_ev::Float64,
    omega_x::Vector{Float64}, 
    omega_ir::Vector{Float64},
    mat::Matrix{Float64},
    outfile::String,
    color_range::Tuple{Float64, Float64},
)
##================================================================##
    gamma = 0.5 
    nx, ny = size(mat)

    x = omega_x  .* au_to_ev
    y = omega_ir .* au_to_ev

    dx = x[2] - x[1]
    dy = y[2] - y[1]

    xmin, xmax = minimum(x) - dx/2, maximum(x) + dx/2
    ymin, ymax = minimum(y) - dy/2, maximum(y) + dy/2

    # Original probability range
    pmin, pmax = color_range

    # Transform matrix and color range using gamma scaling.
    # This stretches low/mid probabilities and compresses high probabilities.
    mat_plot = mat .^ gamma
    color_range_plot = (pmin^gamma, pmax^gamma)

    f = Figure(size = (700, 700), backgroundcolor = :white)

    ax = Axis(
        f[1, 1], 
        xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", 
        ylabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}",
        xlabelsize = 26,
        ylabelsize = 26, 
        xlabelpadding = 8,
        ylabelpadding = 8,
        xticklabelsize = 18,
        yticklabelsize = 18,
        xticksvisible = true,
        yticksvisible = true,
        xticklabelsvisible = true,
        yticklabelsvisible = true,
        xgridvisible = false,
        ygridvisible = false,
        aspect = AxisAspect(1),
        spinewidth = 2.5
    )

    hm = heatmap!(
        ax,
        x,
        y,
        mat_plot,
        colormap = :matter,
        colorrange = color_range_plot,
        interpolate = false
    )

    tightlimits!(ax)

    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)

    x_mesh = collect(LinRange(xmin, xmax, nx + 1))
    y_mesh = collect(LinRange(ymin, ymax, ny + 1))

    for xv in x_mesh
        lines!(ax, [xv, xv], [ymin, ymax], color = (:brown, 0.55), linewidth = 0.35)
    end

    for yv in y_mesh
        lines!(ax, [xmin, xmax], [yv, yv], color = (:brown, 0.55), linewidth = 0.35)
    end

    # Choose colorbar ticks in original probability units
    p_ticks = collect(range(pmin, pmax, length = 5))

    # Tick positions must be transformed, but labels remain original values
    tick_positions = p_ticks .^ gamma
    tick_labels = [string(round(t, digits = 2)) for t in p_ticks]

    Colorbar(
        f[1, 2],
        hm,
        label = L"P(\omega_{\mathrm{IR}}\mid\omega_{\mathrm{X}})",
        labelsize = 18,
        width = 18,
        height = Relative(0.84),
        valign = :center,
        ticks = (tick_positions, tick_labels),
        ticklabelsvisible = true,
        ticksvisible = true
    )

    save(outfile, f, px_per_unit = 3)

    return nothing
end
##================================================================##
function heatmap_plot_matrix_representation_v3(au_to_ev::Float64, omega_x::Vector{Float64}, 
                                omega_ir::Vector{Float64}, mat::Matrix{Float64}, outfile::String,
                                color_range::Tuple{Float64, Float64})
##================================================================##
    nx, ny = size(mat)

    x = omega_x  .* au_to_ev
    y = omega_ir .* au_to_ev

    dx = x[2] - x[1]
    dy = y[2] - y[1]

    xmin, xmax = minimum(x) - dx/2, maximum(x) + dx/2
    ymin, ymax = minimum(y) - dy/2, maximum(y) + dy/2

    vmin, vmax = color_range
    @assert vmin < 0.0 < vmax "For the interference plot, color_range must include zero."

    zero_pos = (0.0 - vmin) / (vmax - vmin)
    cb_ticks = [vmin, 0.0, 0.2, 0.4, 0.6, 0.8, vmax]
    cb_ticks = unique(round.(cb_ticks, digits = 2))
    cb_ticks = filter(t -> vmin <= t <= vmax, cb_ticks)

    positions = [
                0.0,
                zero_pos * 0.65,
                zero_pos,
                min(zero_pos + 0.12, 1.0),
                min(zero_pos + 0.35, 1.0),
                min(zero_pos + 0.62, 1.0),
                1.0
            ]
    interference_cmap = cgrad(
                            [
                                :steelblue4,
                                :lightskyblue,
                                :white,
                                :wheat,
                                :orange,
                                :deeppink3,
                                :purple4
                            ],
                            positions
                        )

    f = Figure(size = (650, 600),
            backgroundcolor = :white,
            figure_padding = (5, 5, 5, 5)   # left, right, bottom, top
            )
    ax = Axis(f[1, 1], 
            xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", 
            # ylabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}",
            xlabelsize = 28, ylabelsize = 28, 
            xlabelpadding = 5, 
            # ylabelpadding = 4,
            xticklabelsize = 18, yticklabelsize = 18,
            xticksvisible = true, yticksvisible = true,
            xticklabelsvisible = true, yticklabelsvisible = true,
            xgridvisible = false, ygridvisible = false,
            aspect = AxisAspect(1), spinewidth = 2.5)

    hm = heatmap!(ax, x, y, mat, 
                # colormap = :balance,
                colormap = interference_cmap,
                # colormap = :seismic,
                colorrange = color_range,
                # colormap = :Reds,
                # colormap = cgrad([:white, :linen, :moccasin, :orange, :firebrick, :darkred]), 
                interpolate = false)

    tightlimits!(ax)

    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)

    x_mesh = collect(LinRange(xmin, xmax, nx + 1))
    y_mesh = collect(LinRange(ymin, ymax, ny + 1))

    for xv in x_mesh
        lines!(ax, [xv, xv], [ymin, ymax], color = (:brown, 0.55), linewidth = 0.3)
    end

    for yv in y_mesh
        lines!(ax, [xmin, xmax], [yv, yv], color = (:brown, 0.55), linewidth = 0.3)
    end

    
    Colorbar(f[1, 2], hm,
            # label = L"I(\omega_{\mathrm{IR}}\mid\omega_{\mathrm{X}})",
            # labelsize = 24,
            width = 18,
            height = Relative(0.99),
            valign = :center,
            ticks = cb_ticks,
            ticklabelsvisible = true,
            ticksvisible = true
        )

    colgap!(f.layout, 7)
    rowgap!(f.layout, 5)
    resize_to_layout!(f)
    save(outfile, f, px_per_unit = 3)
    return nothing
end
##================================================================##
function heatmap_plot_matrix_representation_v2(au_to_ev::Float64, omega_x::Vector{Float64}, 
                                omega_ir::Vector{Float64}, mat::Matrix{Float64}, outfile::String,
                                color_range::Tuple{Float64, Float64})
##================================================================##
    nx, ny = size(mat)

    x = omega_x  .* au_to_ev
    y = omega_ir .* au_to_ev

    dx = x[2] - x[1]
    dy = y[2] - y[1]

    xmin, xmax = minimum(x) - dx/2, maximum(x) + dx/2
    ymin, ymax = minimum(y) - dy/2, maximum(y) + dy/2

    f = Figure(size = (650, 600),
            backgroundcolor = :white,
            figure_padding = (5, 5, 5, 5)   # left, right, bottom, top
            )

    ax = Axis(f[1, 1], 
            xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", 
            ylabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}",
            xlabelsize = 28, ylabelsize = 28, 
            xlabelpadding = 5, 
            ylabelpadding = 4,
            xticklabelsize = 18, yticklabelsize = 18,
            xticksvisible = true, yticksvisible = true,
            xticklabelsvisible = true, yticklabelsvisible = true,
            xgridvisible = false, ygridvisible = false,
            aspect = AxisAspect(1), spinewidth = 2.5)

    hm = heatmap!(ax, x, y, mat, 
                # colormap = :heat,
                # colormap = :amp,
                # colormap = :dense,
                colormap = :matter,
                colorrange = color_range,
                # colormap = :Reds,
                # colormap = cgrad([:white, :linen, :moccasin, :orange, :firebrick, :darkred]), 
                interpolate = false)

    tightlimits!(ax)

    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)

    x_mesh = collect(LinRange(xmin, xmax, nx + 1))
    y_mesh = collect(LinRange(ymin, ymax, ny + 1))

    for xv in x_mesh
        lines!(ax, [xv, xv], [ymin, ymax], color = (:brown, 0.55), linewidth = 0.3)
    end

    for yv in y_mesh
        lines!(ax, [xmin, xmax], [yv, yv], color = (:brown, 0.55), linewidth = 0.3)
    end

    
    Colorbar(f[1, 2], hm, 
            # label = L"P(\omega_{\mathrm{IR}}\mid\omega_{\mathrm{X}})",
            # labelsize = 24, 
            width = 18, 
            height = Relative(0.99), 
            # colorrange = (0.0, quantile(vec(mat), 0.99)),
            valign = :center, ticklabelsvisible = true,  ticksvisible = true, 
            # tellheight = true
            )

    colgap!(f.layout, 7)
    rowgap!(f.layout, 5)
    resize_to_layout!(f)

    save(outfile, f, px_per_unit = 3)
    return nothing
end
##================================================================##
function heatmap_plot_matrix_representation_v1(au_to_ev::Float64, omega_x::Vector{Float64}, 
                                omega_ir::Vector{Float64}, mat::Matrix{Float64}, outfile::String)
##================================================================##
    nx, ny = size(mat)

    x = omega_x  .* au_to_ev
    y = omega_ir .* au_to_ev

    dx = x[2] - x[1]
    dy = y[2] - y[1]

    xmin, xmax = minimum(x) - dx/2, maximum(x) + dx/2
    ymin, ymax = minimum(y) - dy/2, maximum(y) + dy/2

    f = Figure(size = (700, 700), backgroundcolor = :white)

    ax = Axis(f[1, 1], 
            xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", 
            ylabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}",
            xlabelsize = 26, ylabelsize = 26, 
            xlabelpadding = 8, ylabelpadding = 8,
            xticklabelsize = 18, yticklabelsize = 18,
            xticksvisible = true, yticksvisible = true,
            xticklabelsvisible = true, yticklabelsvisible = true,
            xgridvisible = false, ygridvisible = false,
            aspect = AxisAspect(1), spinewidth = 2.5)

    hm = heatmap!(ax, x, y, mat, colormap = :matter, interpolate = false)

    tightlimits!(ax)

    xlims!(ax, xmin, xmax)
    ylims!(ax, ymin, ymax)

    x_mesh = collect(LinRange(xmin, xmax, nx + 1))
    y_mesh = collect(LinRange(ymin, ymax, ny + 1))

    for xv in x_mesh
        lines!(ax, [xv, xv], [ymin, ymax], color = (:brown, 0.55), linewidth = 0.35)
    end

    for yv in y_mesh
        lines!(ax, [xmin, xmax], [yv, yv], color = (:brown, 0.55), linewidth = 0.35)
    end

    Colorbar(f[1, 2], hm, label = L"P(\omega_{\mathrm{X}},\,\omega_{\mathrm{IR}})",
            labelsize = 18, width = 18, 
            height = Relative(0.84), 
            colorrange = (0.0, quantile(vec(mat), 0.99)),
            valign = :center, ticklabelsvisible = true,  ticksvisible = true, 
            tellheight = true)

    save(outfile, f, px_per_unit = 3)
    return nothing
end
##================================================================##
function marginal_ir_helper(mat)
##================================================================##
    nx, ny  = size(mat)
    marg_ir = zeros(ny)

    for ir=1:ny
        marg_ir[ir] = sum(mat[:, ir])  
    end

    return marg_ir
end
##================================================================##
function get_marginal_ir_binned_eem(matrix_c, matrix_q)
##================================================================##
    num_x_c, num_ir_c = size(matrix_c)
    num_x_q, num_ir_q = size(matrix_q)

    marg_ir_c = zeros(num_ir_c)
    marg_ir_q = zeros(num_ir_q)
    @assert size(matrix_c, 2) == size(matrix_q, 2) "The next for-loop is only for equal number of columns in matrix_c and matrix_q"
    
    for ir=1:num_ir_c
        marg_ir_c[ir] = sum(matrix_c[:, ir])  
        marg_ir_q[ir] = sum(matrix_q[:, ir])
    end

    return marg_ir_c, marg_ir_q
end
##================================================================##
function marginal_ir_binned_eem(au_to_ev, omega_x_c, omega_ir_c, matrix_c, omega_x_q, omega_ir_q, matrix_q, plt_path)
##================================================================##
    num_x_c, num_ir_c = size(matrix_c)
    num_x_q, num_ir_q = size(matrix_q)

    marg_ir_c = zeros(num_ir_c)
    marg_ir_q = zeros(num_ir_q)
    @assert size(matrix_c, 2) == size(matrix_q, 2) "The next for-loop is only for equal number of columns in matrix_c and matrix_q"
    
    for ir=1:num_ir_c
        marg_ir_c[ir] = sum(matrix_c[:, ir])  
        marg_ir_q[ir] = sum(matrix_q[:, ir])
    end

    @show sum(marg_ir_c), sum(marg_ir_q)
    omega_ir_c = omega_ir_c .* au_to_ev
    omega_ir_q = omega_ir_q .* au_to_ev

    f  = Figure(size = (1200, 900))
    ax = Axis(f[1,1], 
        xlabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}", ylabel = L"M_{\mathrm{IR}}(\omega_{\mathrm{IR}})",
        xlabelsize = 36, ylabelsize = 36,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = false, yminorgridvisible = false,
        xticksvisible = true, yticksvisible = true,
        xticklabelsize = 24, yticklabelsize = 24,
        xminorticksvisible = true, yminorticksvisible = true,
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 0.66, ytickalign = 0.66,
        xminortickalign = 1, yminortickalign = 1,
        xticksize = 15, yticksize = 15,
        xminorticksize = 5, yminorticksize = 5, 
        xlabelpadding = 18, ylabelpadding = 18,
        )

    lines!(ax, omega_ir_c, marg_ir_c, linewidth = 3.0, color = :darkred, label = "Incoherent")
    lines!(ax, omega_ir_q, marg_ir_q, linewidth = 2.8, color = :darkblue, label = "Coherent")

    leg = Legend(f[1,1], ax, patchsize = (30,20), padding = (15, 15, 15, 15), margin = (15,15,15,15), 
        halign = :right, valign = :top, tellwidth = false,
        labelsize = 28, markersize = 26)

    ax_inset = Axis(f[1, 1],
                width=Relative(0.2),
                height=Relative(0.2),
                halign=0.1,
                valign=0.9,
                # halign=0.9,
                # valign=0.75,
                xticklabelsize = 14,
                yticklabelsize = 14,
                )

    line_inset = lines!(ax_inset, omega_ir_c, marg_ir_c, linewidth = 2.0, color=:darkred)
    translate!(ax_inset.blockscene, 0, 0, 150)
    
    save(plt_path, f, px_per_unit = 2)
    println("marginal_plt: ", plt_path)

    return nothing
end
##================================================================##
function surface_3d_plot_conditional_binned_eem(au_to_ev::Float64, x::Vector{Float64}, y::Vector{Float64}, z::Matrix{Float64}, outfile::String)
##================================================================##
    nx, ny = size(z)
    # @assert length(x) == nx "Dimension mismatch in x length!"
    # @assert length(y) == ny "Dimension mismatch in y length!"

    x = x .* au_to_ev
    y = y .* au_to_ev

    f  = Figure(size = (1200, 900))
    ax = Axis3(f[1, 1], 
                xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}",  
                ylabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}",
                zlabel = L"P(\omega_{\mathrm{IR}}\mid\omega_{\mathrm{X}})",
                xlabelsize = 36, ylabelsize = 36, zlabelsize = 36,
                xticksvisible = true, yticksvisible = true, zticksvisible = true,
                xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22,
                xgridvisible = true, ygridvisible = true, zgridvisible = true, 
                protrusions = (30, 0, 20, 0),  # left, right, bottom, top
                xlabeloffset = 55, ylabeloffset = 55, zlabeloffset = 60,
                # azimuth = 45 * pi / 180,
                elevation = 25 * pi / 180,
                )

    # xlim = 45:75
    # ylim = 1:ny
    # x = x[xlim]
    # y = y[ylim]
    # z = z[xlim,ylim]
    sf = surface!(ax, x, y, z, colormap = :YlOrRd, transparency = true, 
                alpha = 0.75,
                )  #matter #OrRd #YlOrRd
  
    save(outfile, f, px_per_unit = 3)
    return nothing
end
##================================================================##
function max_contributer_v2(omega_x::Vector{Float64}, omega_ir::Vector{Float64}, 
                            matrix::Matrix{Float64}, outfile_1::String, outfile_2::String)
##================================================================##
    max_val, pair = findmax(matrix)
    ix = pair[1]
    ir = pair[2]

    omega_x  = omega_x .* au_to_ev
    omega_ir = omega_ir .* au_to_ev
    val_x    = round(Int, omega_x[ix])
    val_ir   = round(omega_ir[ir]; digits=4)


    fig_1 = Figure(size = (1200, 900))
    ax1   = Axis(fig_1[1,1], 
            xlabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}", 
            ylabel = LaTeXString("P(\\omega_{\\mathrm{X}}=$(val_x) \\mathrm{eV},\\,\\omega_{\\mathrm{IR}})"),
            xlabelsize = 36, ylabelsize = 36,
            xgridvisible = true, ygridvisible = true,
            xminorgridvisible = false, yminorgridvisible = false,
            xticksvisible = true, yticksvisible = true,
            xticklabelsize = 24, yticklabelsize = 24,
            xminorticksvisible = true, yminorticksvisible = true,
            xticksmirrored = true, yticksmirrored = true,
            xtickalign = 0.66, ytickalign = 0.66,
            xminortickalign = 1, yminortickalign = 1,
            xticksize = 15, yticksize = 15,
            xminorticksize = 5, yminorticksize = 5, 
            xlabelpadding = 18, ylabelpadding = 18,
            )
    

    lines!(ax1, omega_ir, matrix[ix, :], linewidth=3.0, color = :darkred)
    save(outfile_1, fig_1, px_per_unit = 2)

    # for i=1:length(omega_ir)
    #     println("$(round(omega_ir[i]; digits=4)), $(round(matrix[ix, i]; digits=4))")
    # end

    #-------
    fig_2 = Figure(size = (1200, 900))
    ax2   = Axis(fig_2[1,1], 
        xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}", 
        ylabel = LaTeXString("P(\\omega_{\\mathrm{X}},\\,\\omega_{\\mathrm{IR}}=$(val_ir) \\mathrm{eV})"),
        xlabelsize = 36, ylabelsize = 36,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = false, yminorgridvisible = false,
        xticksvisible = true, yticksvisible = true,
        xticklabelsize = 24, yticklabelsize = 24,
        xminorticksvisible = true, yminorticksvisible = true,
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 0.66, ytickalign = 0.66,
        xminortickalign = 1, yminortickalign = 1,
        xticksize = 15, yticksize = 15,
        xminorticksize = 5, yminorticksize = 5, 
        xlabelpadding = 18, ylabelpadding = 18,
        )

    xlim = 45:75
    omega_x = omega_x[xlim]
    matrix  = matrix[:, ir][xlim]
    # matrix = matrix[:, ir+1]
    lines!(ax2, omega_x, matrix, linewidth=3.0, color = :darkred)

    save(outfile_2, fig_2, px_per_unit = 2)

    # for i=1:length(omega_x)
    #     # println("$(round(omega_x[i])), $(round(matrix[i]; digits=4))")
    #     println("$(round(omega_x[i])), $(matrix[i])")
    # end
    return nothing
end
##================================================================##
function max_contributer(au_to_ev::Float64, omega_x::Vector{Float64}, omega_ir::Vector{Float64}, matrix::Matrix{Float64}, outfile::String)
##================================================================##
    num_x, num_ir = size(matrix)
    max_val, pair = findmax(matrix)
    ix = pair[1]
    ir = pair[2]

    omega_x  = omega_x .* au_to_ev
    omega_ir = omega_ir .* au_to_ev
    val_x    = round(Int, omega_x[ix])
    val_ir   = round(omega_ir[ir]; digits=4)


    fig  = Figure(size = (1200, 900))
    ax1  = Axis(fig[1,1], 
            xlabel = L"\omega_{IR}\ \mathrm{(eV)}", 
            ylabel = LaTeXString("P(\\omega_X=$(val_x) \\mathrm{(eV)},\\,\\omega_{IR})"),
            xlabelsize = 26, ylabelsize = 26,
            xgridvisible = true, ygridvisible = true,
            xminorgridvisible = false, yminorgridvisible = false,
            xticksvisible = true, yticksvisible = true,
            xticklabelsize = 20, yticklabelsize = 20,
            xminorticksvisible = true, yminorticksvisible = true,
            xticksmirrored = true, yticksmirrored = true,
            xtickalign = 0.66, ytickalign = 0.66,
            xminortickalign = 1, yminortickalign = 1,
            xticksize = 15, yticksize = 15,
            xminorticksize = 5, yminorticksize = 5, 
            ylabelpadding = 15,
            )
    

    lines!(ax1, omega_ir, matrix[ix, :], linewidth=2.5, color = :darkred)

    #-------
    ax2  = Axis(fig[2,1], 
        xlabel = L"\omega_X\ \mathrm{(eV)}", 
        ylabel = LaTeXString("P(\\omega_X,\\,\\omega_{IR}=$(val_ir) \\mathrm{(eV)})"),
        xlabelsize = 26, ylabelsize = 26,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = false, yminorgridvisible = false,
        xticksvisible = true, yticksvisible = true,
        xticklabelsize = 20, yticklabelsize = 20,
        xminorticksvisible = true, yminorticksvisible = true,
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 0.66, ytickalign = 0.66,
        xminortickalign = 1, yminortickalign = 1,
        xticksize = 15, yticksize = 15,
        xminorticksize = 5, yminorticksize = 5,
        ylabelpadding = 15, 
        )

    xlim = 45:75
    omega_x = omega_x[xlim]
    matrix  = matrix[:, ir][xlim]
    lines!(ax2, omega_x, matrix, linewidth=2.5, color = :darkred)

    save(outfile, fig)

    return nothing
end
##================================================================##
function surface_3d_plot_binned_eem(au_to_ev::Float64, x::Vector{Float64}, y::Vector{Float64}, z::Matrix{Float64}, outfile::String)
##================================================================##
    nx, ny = size(z)
    # @assert length(x) == nx "Dimension mismatch in x length!"
    # @assert length(y) == ny "Dimension mismatch in y length!"

    x = x .* au_to_ev
    y = y .* au_to_ev

    f = Figure(size = (1200, 900))
    ax = Axis3(f[1, 1], 
                xlabel = L"\omega_{\mathrm{X}}\ \mathrm{(eV)}",  
                ylabel = L"\omega_{\mathrm{IR}}\ \mathrm{(eV)}",
                zlabel = L"P(\omega_{\mathrm{X}},\,\omega_{\mathrm{IR}})",
                xlabelsize = 36, ylabelsize = 36, zlabelsize = 36,
                xticksvisible = true, yticksvisible = true, zticksvisible = true,
                xticklabelsize = 22, yticklabelsize = 22, zticklabelsize = 22,
                xgridvisible = true, ygridvisible = true, zgridvisible = true, 
                protrusions = (35, 0, 20, 0),  # left, right, bottom, top
                xlabeloffset = 55, ylabeloffset = 55, zlabeloffset = 80,
                # azimuth = -45 * pi / 180,
                elevation = 25 * pi / 180,
                )

    xlim = 45:75
    ylim = 1:ny
    x = x[xlim]
    y = y[ylim]
    z = z[xlim,ylim]
    sf = surface!(ax, x, y, z, colormap = :YlOrRd, transparency = true)  #matter #OrRd #YlOrRd
  
    save(outfile, f, px_per_unit = 3)

    return nothing
end
##================================================================##
function binning_2d(nbins_x::Int, nbins_y::Int, x_vec::Vector{Float64}, y_vec::Vector{Float64}, mat::Matrix{Float64})
##================================================================##
    nx, ny = size(mat)
    x_min, x_max = minimum(x_vec), maximum(x_vec)
    y_min, y_max = minimum(y_vec), maximum(y_vec)
    xbins_range  = collect(LinRange(x_min, x_max, nbins_x+1))
    ybins_range  = collect(LinRange(y_min, y_max, nbins_y+1))

    dx = x_vec[2] - x_vec[1]
    dy = y_vec[2] - y_vec[1]

    # @show dx, dy

    binned_mat   = zeros(nbins_x, nbins_y)
    for ix=1:nbins_x
        x1 = xbins_range[ix]
        x2 = xbins_range[ix+1]
        x_idx = findall(x -> (x >= x1) && (x < x2), x_vec)
        if ix == nbins_x
            x_idx = findall(x -> (x >= x1) && (x <= x2), x_vec)
        end
        for iy=1:nbins_y
            y1 = ybins_range[iy]
            y2 = ybins_range[iy+1]
            y_idx = findall(y -> (y >= y1) && (y < y2), y_vec)
            if iy == nbins_y
                y_idx = findall(y -> (y >= y1) && (y <= y2), y_vec)
            end
            for i in x_idx
                for j in y_idx
                    binned_mat[ix, iy] += mat[i, j] * dx * dy
                end
            end
        end
    end

    # println("Summing up the binned matrix: $(sum(binned_mat)), δω_X = $((xbins_range[2] - xbins_range[1])/2), δω_IR = $((ybins_range[2] - ybins_range[1])/2)") 

    x_centers = (xbins_range[1:end-1] .+ xbins_range[2:end]) ./ 2
    y_centers = (ybins_range[1:end-1] .+ ybins_range[2:end]) ./ 2

    # xlim = 45:75
    # @show xbins_range[44]
    # @show xbins_range[76]
    # @show xbins_range[44] * 27.211
    # @show xbins_range[76] * 27.211

    return x_centers, y_centers, binned_mat
end
##================================================================##
function write_n_max_intensity_to_outfile_v2(au_to_ev::Float64, relative::Float64, x::Vector{Float64}, y::Vector{Float64}, mat::Matrix{Float64}, n::Int, outfile::String)
##================================================================##
    nx      = length(x)
    max_val = maximum(mat)
    fout    = open(outfile, "w")
        title_line = "ω_X (eV) \t"
        for in=1:n
            title_line *= "Peak_$(in) \t ω_IR_$(in) \t"
        end
        write(fout, title_line * "\n")
        for ix=1:nx
            sorted_indices = sortperm(mat[ix, :], rev=true)
            if mat[ix, sorted_indices[1]] < (relative * max_val)
                continue
            else 
                line = "$(x[ix] * au_to_ev)"
                for in=1:n
                    y_idx          = sorted_indices[in]
                    max_intensity  = mat[ix, y_idx]
                    if max_intensity >= (relative * max_val)
                        line *= "\t $max_intensity \t $(y[y_idx] * au_to_ev)"
                    end
                end
                write(fout, line * "\n")
            end
        end
    close(fout)
    return nothing
end
##================================================================##
function write_n_max_intensity_to_outfile(au_to_ev::Float64, x::Vector{Float64}, y::Vector{Float64}, mat::Matrix{Float64}, n::Int, outfile::String)
##================================================================##
    nx = length(x)
    fout = open(outfile, "w")
        title_line = "ω_X (eV) \t"
        for in=1:n
            title_line *= "Peak_$(in) \t ω_IR_$(in) \t"
        end
        write(fout, title_line * "\n")
        for ix=1:nx
            line = "$(x[ix] * au_to_ev) \t"
            for in=1:n
                sorted_indices = sortperm(mat[ix, :], rev=true)
                y_idx          = sorted_indices[in]
                max_intensity  = mat[ix, y_idx]
                line *= "$max_intensity \t $(y[y_idx] * au_to_ev) \t"
            end
            write(fout, line * "\n")
        end
    close(fout)
    return nothing
end
##================================================================##
function get_conditional_mat_v2(matrix::Matrix{Float64})
##================================================================##
    nx, ny          = size(matrix)
    conditional_mat = zeros(nx, ny)

    for ix=1:nx
        int_x_val              = sum(matrix[ix, :])  
        conditional_mat[ix, :] = matrix[ix, :] ./ int_x_val
    end

    @show sum(conditional_mat)
    return conditional_mat
end
##================================================================##
function get_conditional_mat(omega_ir, probs_matrix)
##================================================================##
    num_x, num_ir               = size(probs_matrix)
    signal_omegair_given_omegax = zeros(num_x, num_ir)

    for ix=1:num_x
        int_omega_x = trapz_1D(omega_ir, probs_matrix[ix, :])  
        for ir=1:num_ir
            signal_omegair_given_omegax[ix, ir] = probs_matrix[ix, ir] / int_omega_x  
        end
    end

    return signal_omegair_given_omegax

end
##================================================================##
function writing_max_peak_table_v5(mat_1::Matrix{Float64}, mat_2::Matrix{Float64}, outfile::String)
##================================================================##
    qd_list    = [L"Cd$_{10}$Te$_{10}$", L"Cd$_{20}$Te$_{20}$", L"Cd$_{30}$Te$_{30}$", L"Cd$_{40}$Te$_{40}$", L"Cd$_{50}$Te$_{50}$"]
    fout       = open(outfile, "w")
    title_line = "QD \t Framework \t M(ω_IR) \t ω_IR \t IR_region \n"
    write(fout, title_line)
    n = size(mat_1, 1)
    for i=1:n 
        incoh_ir_val    = (round(mat_1[i,3]; digits=4))
        coh_ir_val      = (round(mat_2[i,3]; digits=4))
        incoh_ir_region = get_ir_region(incoh_ir_val)
        coh_ir_region   = get_ir_region(coh_ir_val)
        line_1 = "$(qd_list[i]) & Incoherent & $(round(mat_1[i,1]; digits=4)) & $incoh_ir_val & $incoh_ir_region \\\\ \n"
        line_2 = "$(qd_list[i]) & Coherent & $(round(mat_2[i,1]; digits=4)) & $coh_ir_val & $coh_ir_region \\\\  \n"
        write(fout, line_1)
        write(fout, line_2)
        write(fout, "\\midrule \n")
    end
    close(fout)
    println("outfile: ", outfile)
    return nothing
end
##================================================================##
function writing_max_peak_table_v4(mat_1::Matrix{Float64}, mat_2::Matrix{Float64}, outfile::String)
##================================================================##
    qd_list    = [L"Cd$_{10}$Te$_{10}$", L"Cd$_{20}$Te$_{20}$", L"Cd$_{30}$Te$_{30}$", L"Cd$_{40}$Te$_{40}$", L"Cd$_{50}$Te$_{50}$"]
    fout       = open(outfile, "w")
    title_line = "QD \t Framework \t P(ω_IR|ω_X) \t ω_X \t ω_IR \t IR_region \n"
    write(fout, title_line)
    n = size(mat_1, 1)
    for i=1:n 
        incoh_ir_val    = (round(mat_1[i,3]; digits=4))
        coh_ir_val      = (round(mat_2[i,3]; digits=4))
        incoh_ir_region = get_ir_region(incoh_ir_val)
        coh_ir_region   = get_ir_region(coh_ir_val)
        line_1 = "$(qd_list[i]) & Incoherent & $(round(mat_1[i,1]; digits=4)) & $(round(Int, mat_1[i,2])) & $incoh_ir_val & $incoh_ir_region \\\\ \n"
        line_2 = "$(qd_list[i]) & Coherent   & $(round(mat_2[i,1]; digits=4)) & $(round(Int, mat_2[i,2])) & $coh_ir_val & $coh_ir_region \\\\ \n"
        write(fout, line_1)
        write(fout, line_2)
        write(fout, "\\midrule \n")
    end
    close(fout)
    println("outfile: ", outfile)
    return nothing
end
##================================================================##
function writing_max_peak_table_v3(mat_1::Matrix{Float64}, mat_2::Matrix{Float64}, outfile::String)
##================================================================##
    qd_list    = [L"Cd$_{10}$Te$_{10}$", L"Cd$_{20}$Te$_{20}$", L"Cd$_{30}$Te$_{30}$", L"Cd$_{40}$Te$_{40}$", L"Cd$_{50}$Te$_{50}$"]
    fout       = open(outfile, "w")
    title_line = "QD \t Formulation \t P(ω_X, ω_IR) \t ω_X \t ω_IR \t IR_region\n"
    write(fout, title_line)
    n = size(mat_1, 1)
    for i=1:n 
        incoh_ir_val    = (round(mat_1[i,3]; digits=4))
        coh_ir_val      = (round(mat_2[i,3]; digits=4))
        incoh_ir_region = get_ir_region(incoh_ir_val)
        coh_ir_region   = get_ir_region(coh_ir_val)
        line_1 = "$(qd_list[i]) & Incoherent & $(round(mat_1[i,1]; digits=4)) & $(round(Int, mat_1[i,2])) & $incoh_ir_val & $incoh_ir_region \\\\ \n"
        line_2 = "$(qd_list[i]) & Coherent  & $(round(mat_2[i,1]; digits=4)) & $(round(Int, mat_2[i,2])) & $coh_ir_val & $coh_ir_region \\\\ \n"
        write(fout, line_1)
        write(fout, line_2)
        write(fout, "\\midrule \n")
    end
    close(fout)
    println("outfile: ", outfile)
    return nothing
end
##================================================================##
function writing_max_peak_table_v2(mat_1::Matrix{Float64}, mat_2::Matrix{Float64}, mat_3::Matrix{Float64}, 
                                    mat_4::Matrix{Float64}, outfile::String)
##================================================================##
    qd_list    = [L"Cd_{10}Te_{10}", L"Cd_{20}Te_{20}", L"Cd_{30}Te_{30}", L"Cd_{40}Te_{40}", L"Cd_{50}Te_{50}"]
    fout       = open(outfile, "w")
    title_line = "QD \t Framework \t Method \t P(ω_IR|ω_X) \t ω_X \t ω_IR \n"
    write(fout, title_line)
    n = size(mat_1, 1)
    for i=1:n 
        line_1 = "$(qd_list[i]) \t Incoherent \t HF \t $(mat_1[i,1]) \t $(mat_1[i,2]) \t $(mat_1[i,3]) \n"
        line_2 = "$(qd_list[i]) \t Incoherent \t DFT \t $(mat_2[i,1]) \t $(mat_2[i,2]) \t $(mat_2[i,3]) \n"
        line_3 = "$(qd_list[i]) \t Coherent \t HF \t $(mat_3[i,1]) \t $(mat_3[i,2]) \t $(mat_3[i,3]) \n"
        line_4 = "$(qd_list[i]) \t Coherent \t DFT \t $(mat_4[i,1]) \t $(mat_4[i,2]) \t $(mat_4[i,3]) \n"
        write(fout, line_1)
        write(fout, line_2)
        write(fout, line_3)
        write(fout, line_4)
    end
    close(fout)
    println("outfile: ", outfile)
    return nothing
end
##================================================================##
function writing_max_peak_table(mat::Matrix{Float64}, outfile::String)
##================================================================##
    qd_list    = [L"Cd_{10}Te_{10}", L"Cd_{20}Te_{20}", L"Cd_{30}Te_{30}", L"Cd_{40}Te_{40}", L"Cd_{50}Te_{50}"]
    fout       = open(outfile, "w")
    title_line = "QD \t P(ω_X, ω_IR) \t ω_X \t ω_IR \n"
    write(fout, title_line)
    n = size(mat, 1)
    for i=1:n 
        line = "$(qd_list[i]) \t $(mat[i,1]) \t $(mat[i,2]) \t $(mat[i,3]) \n"
        write(fout, line)
    end
    close(fout)
    println("outfile: ", outfile)
    return nothing
end
##================================================================##
function single_max_peak_location_v2(x::Vector{Float64}, y::Vector{Float64})
##================================================================##
    au_to_ev = 27.211
    max_peak, loc = findmax(y)
    x_val = x[loc] * au_to_ev 
    return [max_peak, 0.0, x_val]
end
##================================================================##
function single_max_peak_location(x::Vector{Float64}, y::Vector{Float64}, z::Matrix{Float64})
##================================================================##
    max_peak, pair = findmax(z)
    x_val = x[pair[1]] * au_to_ev 
    y_val = y[pair[2]] * au_to_ev
    return [max_peak, x_val, y_val]
end
##================================================================##
function writing_marginal_ir_all_size(x::Vector{Float64}, mat::Matrix{Float64}, outfile::String)
##================================================================##
    au_to_ev = 27.211
    x = x .* au_to_ev

    fout = open(outfile, "w")
    title_line = "omega_IR(eV) \t CdTe_20 \t CdTe_40 \t CdTe_60 \t CdTe_80 \t CdTe_100 \n"
    write(fout, title_line)
    n = length(x)
    for i=1:n
        line = "$(x[i]) \t $(mat[1,i]) \t $(mat[2,i]) \t $(mat[3,i]) \t $(mat[4,i]) \t $(mat[5,i]) \n"
        write(fout, line)
    end
    close(fout)
    println("outfile: ", outfile)
    return nothing
end
##================================================================##
function marginal_x_all_size(omega_x, mat, plt_path)
##================================================================##
    au_to_ev = 27.211
    omega_x = omega_x .* au_to_ev

    f = Figure(size = (1200, 900))
    ax = Axis(f[1,1], #title = "Marginal Probability Distribution",
        xlabel = L"\omega_X\ (\mathrm{eV})", ylabel = L"M_X\ (\omega_X)",
        xlabelsize = 24, ylabelsize = 24,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = false, yminorgridvisible = false,
        xticksvisible = true, yticksvisible = true,
        xticklabelsize = 20, yticklabelsize = 20,
        xminorticksvisible = true, yminorticksvisible = true,
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 0.66, ytickalign = 0.66,
        xminortickalign = 1, yminortickalign = 1,
        xticksize = 15, yticksize = 15,
        xminorticksize = 5, yminorticksize = 5, 
        )
    
    lines!(ax, omega_x, mat[1,:], label = L"Cd_{10}\,Te_{10}")
    lines!(ax, omega_x, mat[2,:], label = L"Cd_{20}Te_{20}")
    lines!(ax, omega_x, mat[3,:], label = L"Cd_{30}Te_{30}")
    lines!(ax, omega_x, mat[4,:], label = L"Cd_{40}Te_{40}")
    lines!(ax, omega_x, mat[5,:], label = L"Cd_{50}Te_{50}")

    leg = Legend(f[1,1], ax, patchsize = (20,20), padding = (10,10,10,10), margin = (10,10,10,10), 
        halign = :right, valign = :top, tellwidth = false,
        labelsize = 20, markersize = 20)    

    save(plt_path, f)
    println("marginal_plt: ", plt_path)

    return nothing
end
##================================================================##
function marginal_ir_all_size(omega_ir, mat, plt_path)
##================================================================##
    au_to_ev = 27.211
    omega_ir = omega_ir .* au_to_ev

    f  = Figure(size = (1200, 900))
    ax = Axis(f[1,1], 
        xlabel = L"\omega_{IR}\ \mathrm{(eV)}", ylabel = L"M_{IR}(\omega_{IR})",
        xlabelsize = 24, ylabelsize = 24,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = false, yminorgridvisible = false,
        xticksvisible = true, yticksvisible = true,
        xticklabelsize = 20, yticklabelsize = 20,
        xminorticksvisible = true, yminorticksvisible = true,
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 0.66, ytickalign = 0.66,
        xminortickalign = 1, yminortickalign = 1,
        xticksize = 15, yticksize = 15,
        xminorticksize = 5, yminorticksize = 5, 
        )
    
    lines!(ax, omega_ir, mat[1,:], label = L"Cd_{10}Te_{10}")
    lines!(ax, omega_ir, mat[2,:], label = L"Cd_{20}Te_{20}")
    lines!(ax, omega_ir, mat[3,:], label = L"Cd_{30}Te_{30}")
    lines!(ax, omega_ir, mat[4,:], label = L"Cd_{40}Te_{40}")
    lines!(ax, omega_ir, mat[5,:], label = L"Cd_{50}Te_{50}")

    leg = Legend(f[1,1], ax, patchsize = (20,20), padding = (10,10,10,10), margin = (10,10,10,10), 
        halign = :right, valign = :top, tellwidth = false,
        labelsize = 20, markersize = 20)    

    save(plt_path, f)
    println("marginal_plt: ", plt_path)

    return nothing
end
##================================================================##
function marginal_ir_pair_hf_and_dft(omega_x1, omega_ir1, mat1, omega_x2, omega_ir2, mat2, plt_path)
##================================================================##
    num_x1, num_ir1 = size(mat1)
    num_x2, num_ir2 = size(mat2)
    marg_ir1 = zeros(num_ir1)
    marg_ir2 = zeros(num_ir2)
    @assert num_ir1 == num_ir2 "The next for-loop is only for equal number of num_ir1 and num_ir2"

    for ir=1:num_ir1
        marg_ir1[ir] = trapz_1D(omega_x1, mat1[:, ir])  # integrate vs ω_X (rows) for each fixed ω_IR (column)
        marg_ir2[ir] = trapz_1D(omega_x2, mat2[:, ir])
    end

    auc = trapz_1D(omega_ir1, marg_ir1)
    println("AUC of marginal_IR1 = $auc")

    au_to_ev = 27.211
    omega_ir1 = omega_ir1 .* au_to_ev
    omega_ir2 = omega_ir2 .* au_to_ev

    f = Figure(size = (1200, 900))
    ax = Axis(f[1,1], #title = "Marginal Probability Distribution",
        xlabel = "Emission Frequency (eV)", ylabel = "IR Marginal Probability Density",
        xlabelsize = 20, ylabelsize = 20,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = false, yminorgridvisible = false,
        xticksvisible = true, yticksvisible = true,
        xticklabelsize = 20, yticklabelsize = 20,
        xminorticksvisible = true, yminorticksvisible = true,
        xticksmirrored = true, yticksmirrored = true,
        xtickalign = 0.66, ytickalign = 0.66,
        xminortickalign = 1, yminortickalign = 1,
        xticksize = 15, yticksize = 15,
        xminorticksize = 5, yminorticksize = 5, 
        )
    
    lines!(ax, omega_ir1, marg_ir1, color = :blue, label = "HF")
    lines!(ax, omega_ir2, marg_ir2, color = :red, label = "DFT")

    leg = Legend(f[1,1], ax, patchsize = (20,20), padding = (10,10,10,10), margin = (10,10,10,10), 
        halign = :right, valign = :top, tellwidth = false,
        labelsize = 20, markersize = 20)

    # ax_inset = Axis(f[1, 1],
    #             width=Relative(0.2),
    #             height=Relative(0.2),
    #             halign=0.1,
    #             valign=0.9,
    #             # title="incoherent"
    #             )

    # line_inset = lines!(ax_inset, omega_ir_c, marg_ir_c, color=:blue)
    # translate!(ax_inset.blockscene, 0, 0, 150)
    

    save(plt_path, f)
    println("marginal_plt: ", plt_path)
    return nothing
end
##================================================================##
function writing_max_peaks_location_to_outfile(n::Int, au_to_ev::Float64, x::Vector{Float64}, y::Vector{Float64}, mat::Matrix{Float64}, outfile::String)
##================================================================##
    max_peak_pairs = get_location_of_peaks(n, mat)
    fout = open(outfile, "w")
        write(fout, "P(ω_X, ω_IR) \t ω_X \t ω_IR \n")
        for i=1:n
            pair = max_peak_pairs[i]
            prob = mat[pair[1],pair[2]]
            line = "$prob \t $(x[pair[1]]*au_to_ev) \t $(y[pair[2]]*au_to_ev) \n"
            write(fout, line)
        end
    close(fout)
    println("outfile: ", outfile)
    return nothing
end
##================================================================##
function surface_3d_plot(x::Vector{Float64}, y::Vector{Float64}, z::Matrix{Float64}, outfile::String)
##================================================================##
    nx, ny = size(z)
    # @assert length(x) == nx "Dimension mismatch in x length!"
    # @assert length(y) == ny "Dimension mismatch in y length!"

    # Creates the figure and sets the size
    f = Figure(size = (1200, 900))

    # Creates the 3D axis with gridlines for publication quality
    ax = Axis3(f[1, 1], 
                xlabel = "Excitation Frequency (eV)", ylabel = "Emission Frequency (eV)", zlabel = "Relative Probability of Emission",
                xlabelsize = 20, ylabelsize = 20, zlabelsize = 20,
                xticksvisible = true, yticksvisible = true, zticksvisible = true,
                xticklabelsize = 15, yticklabelsize = 15, zticklabelsize = 15,
                xgridvisible = true, ygridvisible = true, zgridvisible = true, 
                # ztickformat = "{:.1e}",
                zlabeloffset= 80, protrusions = (50, 5, 50, 5)  # left, right, bottom, top
                )

    # The colormap keyword determines the color gradient for the surface
    sf = surface!(ax, x, y, z, colormap = :YlOrRd)  #matter #OrRd #YlOrRd
  
    # Saves the figure
    save(outfile, f)

    return nothing
end
##================================================================##
function get_location_of_peaks_v2(n::Int, mat::Matrix{Float64})
##================================================================##
    ## This function will return n-pair of (x,y) for n-highest entries mat
    numx, numy = size(mat)
    col_sums   = zeros(numx)
    for iy=1:numy
        col_sums[iy] = sum(mat[:,iy])
    end
    sorted_indices = sortperm(col_sums, rev=true)
    vec = fill((0,0), n)
    for j=1:n 
        idx       = sorted_indices[j]
        val, pair = findmax(mat[:,idx])
        vec[j]    = (idx, pair)
    end
    return vec
    return nothing
end
##================================================================##
function get_location_of_peaks(n::Int, mat::Matrix{Float64})
##================================================================##
    ## This function will return n-pair of (x,y) for n-highest entries mat
    numx, numy = size(mat)
    row_sums   = zeros(numx)
    for ix=1:numx
        row_sums[ix] = sum(mat[ix,:])
    end
    sorted_indices = sortperm(row_sums, rev=true)
    vec = fill((0,0), n)
    for j=1:n 
        idx       = sorted_indices[j]
        val, pair = findmax(mat[idx,:])
        vec[j]    = (idx, pair)
    end
    return vec
end
##================================================================##
function get_n_high_probs(n::Int, mat::Matrix{Float64})
##================================================================##
    ## This function will return n-pair of (x,y) for n-highest entries mat
    numx, numy = size(mat)
    val = zeros(n)
    vec = fill((0,0), n)
    for ix=1:numx
        for iy=1:numy
            entry = mat[ix,iy]
            min_val, idx = findmin(val)
            if entry > min_val
                val[idx] = entry
                vec[idx] = (ix, iy)
            end
        end
    end

    ## To double-check 
    # max_val, idx = findmax(mat)
    # if max_val in val
    #     println("The right max value picked!")
    #     @show max_val
    #     @show idx
    # end
    # @show vec
    return vec
end
##================================================================##
function trapz_1D(x::Vector, f::Vector)
##================================================================##
    num = length(x)
    s   = 0.0e0
    for i=1:num-1
        dx = x[i+1] - x[i]
        s += (0.5 * (f[i+1] + f[i]) * dx)
    end
    return s
end
##================================================================##
function area_under_the_curve_calculator_2D(omega_x, omega_ir, probs_matrix)
##================================================================##
## 2D trapezoidal rule: integrate along ω_IR for each fixed ω_X, then vs ω_X
    num_x   = length(omega_x)
    row_int = zeros(num_x)
    for i=1:num_x
        row_int[i] = trapz_1D(omega_ir, probs_matrix[i,:])
    end
    auc_val = trapz_1D(omega_x, row_int)
    return auc_val
end
##================================================================##
function infile_reader(infile)
##================================================================##
    num        = countlines(infile) - 1  
    omega_x    = zeros(num)
    omega_ir   = zeros(num)
    probs      = zeros(num)
    fin        = open(infile, "r")
    title_line = readline(fin)
    for i=1:num
        cols = split(readline(fin))
        omega_x[i]  = parse(Float64, cols[1])
        omega_ir[i] = parse(Float64, cols[2])
        probs[i]    = parse(Float64, cols[3])
    end

    #-------
    ## Creating EEM using raw data
    unique_omega_x  = sort(unique(omega_x))
    unique_omega_ir = sort(unique(omega_ir))
    num_x  = length(unique_omega_x)
    num_ir = length(unique_omega_ir)
    probs_matrix = zeros(num_x, num_ir)
    counter = 1
    for ix=1:num_x
        for ir=1:num_ir
            probs_matrix[ix,ir] = probs[counter]
            counter += 1
        end
    end

    return num_x, num_ir, unique_omega_x, unique_omega_ir, probs_matrix
end
##================================================================##
function test1A()
##================================================================##
## Binning the EEMs and plotting the unitless EEMs!
    sys_list = [100]
    # gen_dir_path = "./first_order_corr_res/"
    # mkpath(gen_dir_path)

    for (idx, isys) in enumerate(sys_list)
        infiles_dir_path = "./elec_hole_corrected/elec_hole_corrected_raw_data/cdte_$(isys)/"

        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        # infile_incoh_hf  = infiles_dir_path * "cdte$(isys)_hf_classical_prob_res_0.0625.out"
        # infile_coh_hf    = infiles_dir_path * "cdte$(isys)_hf_quantum_prob_res_0.0625.out"
        infile_incoh_dft = infiles_dir_path * "cdte$(isys)_b3lyp_classical_prob_res_0.0625.out"
        infile_coh_dft   = infiles_dir_path * "cdte$(isys)_b3lyp_quantum_prob_res_0.0625.out"

        # num_x_incoh_hf, num_ir_incoh_hf, omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf      = infile_reader(infile_incoh_hf)
        # num_x_coh_hf, num_ir_coh_hf, omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf      = infile_reader(infile_coh_hf)
        num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)
        
        #------- EEMs Normalization ------
        # ## incoherent HF and DFT
        # auc_val_cb_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC before normalization: $auc_val_cb_hf")
        # norm_probs_matrix_incoh_hf = probs_matrix_incoh_hf ./ auc_val_cb_hf
        # # auc_val_ca_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC after normalization: $auc_val_ca_hf")

        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC before normalization: $auc_val_cb_dft")
        norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        # auc_val_ca_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC after normalization: $auc_val_ca_dft")

        # ## coherent HF and DFT
        # auc_val_qb_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf)
        # # println("coherent-HF AUC before normalization: $auc_val_qb_hf")
        # norm_probs_matrix_coh_hf = probs_matrix_coh_hf ./ auc_val_qb_hf
        # # auc_val_qa_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # # println("coherent-HF AUC after normalization: $auc_val_qa_hf")

        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        # println("coherent-DFT AUC before normalization: $auc_val_qb_dft")
        norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        # auc_val_qa_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)
        # println("coherent-DFT AUC after normalization: $auc_val_qa_dft")

        ##------ binning the EEMs ------
        nbins_x  = 200
        nbins_ir = 200
        # binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)  
        @show sum(binned_probs_matrix_coh_dft)
        #----- surface plotting the binned EEMs ------
        # println("Plotting the binned EEMs for system size: $isys")
        # plt_dir_path   = "./first_order_corr_res/cdte_$(isys)/surf_plts/"
        # mkpath(plt_dir_path)
        # plt_path_incoh_hf  = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_incoherent_hf.png"
        # plt_path_coh_hf    = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_coherent_hf.png"
        # plt_path_incoh_dft = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_incoherent_b3lyp.png"
        # plt_path_coh_dft   = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_coherent_b3lyp.png"
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf, plt_path_incoh_hf)
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf, plt_path_coh_hf)
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, plt_path_incoh_dft)
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, plt_path_coh_dft)
    
        ##----- Extra data analysis ------
        println("Probability \t omega_X \t omega_IR")
        mat_t = deepcopy(binned_probs_matrix_incoh_dft)
        xlim = 45:75
        mat_t = mat_t[xlim, :]
        x_range = binned_omega_x_coh_dft[xlim]
        for i=1:10
            val, pair = findmax(mat_t)
            # @show pair
            println("$val \t $(x_range[pair[1]]*au_to_ev) \t $(binned_omega_ir_coh_dft[pair[2]]*au_to_ev)")
            mat_t[pair] = 0.0
        end
        # @show length(binned_probs_matrix_coh_dft[xlim,39])
        # @show mean(binned_probs_matrix_coh_dft[xlim,39])
        # @show median(binned_probs_matrix_coh_dft[xlim,39])

        #----- 2D plotting max peaks of the binned EEMs ------
        println("Plotting the binned EEMs for system size: $isys")
        gen_dir_path   = "./first_order_corr_res/cdte_$(isys)/surf_plts/2d_max_plt/"
        mkpath(gen_dir_path)
        plt_path_incoh_dft_x  = gen_dir_path * "cdte_$(isys)_prob_vs_omegax_incoherent_b3lyp.pdf"
        plt_path_incoh_dft_ir = gen_dir_path * "cdte_$(isys)_prob_vs_omegair_incoherent_b3lyp.pdf"
        plt_path_coh_dft_x  = gen_dir_path * "cdte_$(isys)_prob_vs_omegax_coherent_b3lyp.pdf"
        plt_path_coh_dft_ir = gen_dir_path * "cdte_$(isys)_prob_vs_omegair_coherent_b3lyp.pdf"
        max_contributer_v2(binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, plt_path_incoh_dft_ir, plt_path_incoh_dft_x)
        max_contributer_v2(binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, plt_path_coh_dft_ir, plt_path_coh_dft_x)
        
    end
    return nothing
end
##================================================================##
function test1B()
##================================================================##
## Binning the EEMs and plotting the unitless EEMs!
## Max peaks on 2D plots!
    sys_list = [20]
    # gen_dir_path = "./binned_EEMs/"
    # mkpath(gen_dir_path)
    for (idx, isys) in enumerate(sys_list)
        # gen_dir_path = "./cdte_$(isys)/binned_EEMs/"
        # mkpath(gen_dir_path)

        au_to_ev = 27.211
        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        infile_incoh_hf  = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_hf_incoherent_prob_res_0.0625.out"
        infile_coh_hf  = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_hf_coherent_prob_res_0.0625.out"
        infile_incoh_dft = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_b3lyp_incoherent_prob_res_0.0625.out"
        infile_coh_dft = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_b3lyp_coherent_prob_res_0.0625.out"

        num_x_incoh_hf, num_ir_incoh_hf, omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf      = infile_reader(infile_incoh_hf)
        num_x_coh_hf, num_ir_coh_hf, omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf      = infile_reader(infile_coh_hf)
        # num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        # num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)
        
        #------- EEMs Normalization ------
        ## incoherent HF and DFT
        auc_val_cb_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf)
        # println("incoherent-HF AUC before normalization: $auc_val_cb_hf")
        norm_probs_matrix_incoh_hf = probs_matrix_incoh_hf ./ auc_val_cb_hf
        # auc_val_ca_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # println("incoherent-HF AUC after normalization: $auc_val_ca_hf")

        # auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        # # println("incoherent-DFT AUC before normalization: $auc_val_cb_dft")
        # norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        # # auc_val_ca_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # # println("incoherent-DFT AUC after normalization: $auc_val_ca_dft")

        ## coherent HF and DFT
        auc_val_qb_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf)
        # println("coherent-HF AUC before normalization: $auc_val_qb_hf")
        norm_probs_matrix_coh_hf = probs_matrix_coh_hf ./ auc_val_qb_hf
        # auc_val_qa_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # println("coherent-HF AUC after normalization: $auc_val_qa_hf")

        # auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        # # println("coherent-DFT AUC before normalization: $auc_val_qb_dft")
        # norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        # # auc_val_qa_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)
        # # println("coherent-DFT AUC after normalization: $auc_val_qa_dft")

        ##------ binning the EEMs ------
        nbins_x  = 200
        nbins_ir = 200
        binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)  
        
        #----- 2D plotting max peaks of the binned EEMs ------
        println("Plotting the binned EEMs for system size: $isys")
        gen_dir_path   = "./binned_EEMs/cdte_$(isys)/max_peaks_2d_plots/zoomed_in/"
        # mkpath(gen_dir_path)
        plt_path_incoh_hf_x   = gen_dir_path * "cdte_$(isys)_prob_vs_omegax_incoherent_hf.pdf"
        plt_path_incoh_hf_ir  = gen_dir_path * "cdte_$(isys)_prob_vs_omegair_incoherent_hf.pdf"
        plt_path_coh_hf_x   = gen_dir_path * "cdte_$(isys)_prob_vs_omegax_coherent_hf.pdf"
        plt_path_coh_hf_ir  = gen_dir_path * "cdte_$(isys)_prob_vs_omegair_coherent_hf.pdf"
        # plt_path_incoh_dft_x  = gen_dir_path * "cdte_$(isys)_prob_vs_omegax_incoherent_b3lyp.pdf"
        # plt_path_incoh_dft_ir = gen_dir_path * "cdte_$(isys)_prob_vs_omegair_incoherent_b3lyp.pdf"
        # plt_path_coh_dft_x  = gen_dir_path * "cdte_$(isys)_prob_vs_omegax_coherent_b3lyp.pdf"
        # plt_path_coh_dft_ir = gen_dir_path * "cdte_$(isys)_prob_vs_omegair_coherent_b3lyp.pdf"
        # max_contributer_v2(au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf, plt_path_incoh_hf_ir, plt_path_incoh_hf_x)
        # max_contributer_v2(au_to_ev, binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf, plt_path_coh_hf_ir, plt_path_coh_hf_x)
        max_contributer_v2(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, plt_path_incoh_dft_ir, plt_path_incoh_dft_x)
        # max_contributer_v2(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, plt_path_coh_dft_ir, plt_path_coh_dft_x)
        
        #----- printing the max peaks of the binned EEMs ------
        # println("Max peaks for the binned EEMs for system size: $isys")
        # n = 10
        # max_peaks_indices_incoh_hf = sortperm(vec(binned_probs_matrix_incoh_hf), rev=true)[1:n]
        # # max_peaks_indices_coh_hf = sortperm(vec(binned_probs_matrix_coh_hf), rev=true)[1:5]
        # # max_peaks_indices_incoh_dft = sortperm(vec(binned_probs_matrix_incoh_dft), rev=true)[1:5]
        # # max_peaks_indices_coh_dft = sortperm(vec(binned_probs_matrix_coh_dft), rev=true)[1:5]

        # println("----- Incoherent HF max peaks -----")
        # inds = CartesianIndices(binned_probs_matrix_incoh_hf)[max_peaks_indices_incoh_hf]
        # for i=1:n
        #     idx = inds[i]
        #     prob = binned_probs_matrix_incoh_hf[idx]
        #     x_idx, y_idx = Tuple(idx)
        #     println("  Peak $i: (probability = $(binned_probs_matrix_incoh_hf[idx]), ω_X = $(binned_omega_x_incoh_hf[x_idx]*au_to_ev), ω_IR = $(binned_omega_ir_incoh_hf[y_idx]*au_to_ev))")
        # end
    end
    return nothing
end
##================================================================##
function test1C()
##================================================================##
## Binning the EEMs and plotting the unitless EEMs!
## cEEM 3D surface plots for the binned EEMs!
    sys_list = [20, 40, 60, 80, 100]

    for (idx, isys) in enumerate(sys_list)
        infiles_dir_path = "./elec_hole_corrected/elec_hole_corrected_raw_data/cdte_$(isys)/"

        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        # infile_incoh_hf  = infiles_dir_path * "cdte$(isys)_hf_classical_prob_res_0.0625.out"
        # infile_coh_hf    = infiles_dir_path * "cdte$(isys)_hf_quantum_prob_res_0.0625.out"
        infile_incoh_dft = infiles_dir_path * "cdte$(isys)_b3lyp_classical_prob_res_0.0625.out"
        infile_coh_dft   = infiles_dir_path * "cdte$(isys)_b3lyp_quantum_prob_res_0.0625.out"

        # num_x_incoh_hf, num_ir_incoh_hf, omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf      = infile_reader(infile_incoh_hf)
        # num_x_coh_hf, num_ir_coh_hf, omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf      = infile_reader(infile_coh_hf)
        num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)
        
        #------- EEMs Normalization ------
        ## incoherent HF and DFT
        # auc_val_cb_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC before normalization: $auc_val_cb_hf")
        # norm_probs_matrix_incoh_hf = probs_matrix_incoh_hf ./ auc_val_cb_hf
        # # auc_val_ca_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC after normalization: $auc_val_ca_hf")

        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC before normalization: $auc_val_cb_dft")
        norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        # auc_val_ca_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC after normalization: $auc_val_ca_dft")

        # # ## coherent HF and DFT
        # auc_val_qb_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf)
        # # println("coherent-HF AUC before normalization: $auc_val_qb_hf")
        # norm_probs_matrix_coh_hf = probs_matrix_coh_hf ./ auc_val_qb_hf
        # # auc_val_qa_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # # println("coherent-HF AUC after normalization: $auc_val_qa_hf")

        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        # println("coherent-DFT AUC before normalization: $auc_val_qb_dft")
        norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        # auc_val_qa_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)
        # println("coherent-DFT AUC after normalization: $auc_val_qa_dft")

        ##------ binning the EEMs ------
        nbins_x  = 200
        nbins_ir = 200
        # binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)  
        
        #----- 3D plotting conditional surfaces of the binned EEMs -------
        println("Plotting the binned EEMs for system size: $isys")
        plt_dir_path   = "./first_order_corr_res/cdte_$(isys)/ceem_plts/"
        # mkpath(plt_dir_path)
        plt_path_incoh_hf  = plt_dir_path * "cdte_$(isys)_3d_conditional_incoherent_hf.png"
        plt_path_coh_hf    = plt_dir_path * "cdte_$(isys)_3d_conditional_coherent_hf.png"
        plt_path_incoh_dft = plt_dir_path * "cdte_$(isys)_3d_conditional_incoherent_b3lyp.png"
        plt_path_coh_dft   = plt_dir_path * "cdte_$(isys)_3d_conditional_coherent_b3lyp.png"
        # conditional_mat_incoh_hf  = get_conditional_mat_v2(binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf)
        # conditional_mat_coh_hf  = get_conditional_mat_v2(binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf)
        conditional_mat_incoh_dft = get_conditional_mat_v2(binned_probs_matrix_incoh_dft)
        conditional_mat_coh_dft = get_conditional_mat_v2(binned_probs_matrix_coh_dft)
        # surface_3d_plot_conditional_binned_eem(au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, conditional_mat_incoh_hf, plt_path_incoh_hf)
        # surface_3d_plot_conditional_binned_eem(au_to_ev, binned_omega_x_coh_hf, binned_omega_ir_coh_hf, conditional_mat_coh_hf, plt_path_coh_hf)
        # surface_3d_plot_conditional_binned_eem(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, conditional_mat_incoh_dft, plt_path_incoh_dft)
        surface_3d_plot_conditional_binned_eem(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, conditional_mat_coh_dft, plt_path_coh_dft)

    end
    return nothing
end
##================================================================##
function test1D()
##================================================================##
## Binning the EEMs and plotting the unitless EEMs!
## Marginals plots for the binned EEMs!
    sys_list = [100]
    
    for (idx, isys) in enumerate(sys_list)
        infiles_dir_path = "./elec_hole_corrected/elec_hole_corrected_raw_data/cdte_$(isys)/"

        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        # infile_incoh_hf  = infiles_dir_path * "cdte$(isys)_hf_classical_prob_res_0.0625.out"
        # infile_coh_hf    = infiles_dir_path * "cdte$(isys)_hf_quantum_prob_res_0.0625.out"
        infile_incoh_dft = infiles_dir_path * "cdte$(isys)_b3lyp_classical_prob_res_0.0625.out"
        infile_coh_dft   = infiles_dir_path * "cdte$(isys)_b3lyp_quantum_prob_res_0.0625.out"

        # num_x_incoh_hf, num_ir_incoh_hf, omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf      = infile_reader(infile_incoh_hf)
        # num_x_coh_hf, num_ir_coh_hf, omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf      = infile_reader(infile_coh_hf)
        num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)

        #------- EEMs Normalization ------
        # ## incoherent HF and DFT
        # auc_val_cb_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC before normalization: $auc_val_cb_hf")
        # norm_probs_matrix_incoh_hf = probs_matrix_incoh_hf ./ auc_val_cb_hf
        # # auc_val_ca_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC after normalization: $auc_val_ca_hf")

        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC before normalization: $auc_val_cb_dft")
        norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        # auc_val_ca_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC after normalization: $auc_val_ca_dft")

        # # ## coherent HF and DFT
        # auc_val_qb_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf)
        # # println("coherent-HF AUC before normalization: $auc_val_qb_hf")
        # norm_probs_matrix_coh_hf = probs_matrix_coh_hf ./ auc_val_qb_hf
        # # auc_val_qa_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # # println("coherent-HF AUC after normalization: $auc_val_qa_hf")

        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        # println("coherent-DFT AUC before normalization: $auc_val_qb_dft")
        norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        # auc_val_qa_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)
        # println("coherent-DFT AUC after normalization: $auc_val_qa_dft")

        ##------ binning the EEMs ------
        nbins_x  = 200
        nbins_ir = 200
        # binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf   = binning_2d(nbins_x, nbins_ir, omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf         = binning_2d(nbins_x, nbins_ir, omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft       = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)  
        
        #----- 2D plotting marginals of the binned EEMs -------
        println("Plotting for system size: $isys")
        plt_dir_path   = "./first_order_corr_res/cdte_$(isys)/marg_plts/"
        # mkpath(plt_dir_path)
        plt_path_hf  = plt_dir_path * "cdte_$(isys)_marginal_ir_hf.pdf"
        plt_path_dft = plt_dir_path * "cdte_$(isys)_marginal_ir_b3lyp.pdf"
         
        # marginal_ir_binned_eem(au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf, 
        #                         binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf, plt_path_hf)
        marginal_ir_binned_eem(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, 
                                binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, plt_path_dft)
    end
    return nothing
end
##================================================================##
function test1E()
##================================================================##
## Binning the EEMs and plotting the unitless EEMs!
## extracting tables data from the binned EEMs!

    sys_list  = [20, 40, 60, 80, 100]
    num_sys   = length(sys_list)
    num_cols  = 3
    mat_incoh_hf  = zeros(num_sys, num_cols)
    mat_incoh_dft = zeros(num_sys, num_cols)
    mat_coh_hf  = zeros(num_sys, num_cols)
    mat_coh_dft = zeros(num_sys, num_cols)

    for (idx, isys) in enumerate(sys_list)

        infiles_dir_path = "./elec_hole_corrected/elec_hole_corrected_raw_data/cdte_$(isys)/"

        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        # infile_incoh_hf  = infiles_dir_path * "cdte$(isys)_hf_classical_prob_res_0.0625.out"
        # infile_coh_hf    = infiles_dir_path * "cdte$(isys)_hf_quantum_prob_res_0.0625.out"
        infile_incoh_dft = infiles_dir_path * "cdte$(isys)_b3lyp_classical_prob_res_0.0625.out"
        infile_coh_dft   = infiles_dir_path * "cdte$(isys)_b3lyp_quantum_prob_res_0.0625.out"

        # num_x_incoh_hf, num_ir_incoh_hf, omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf      = infile_reader(infile_incoh_hf)
        # num_x_coh_hf, num_ir_coh_hf, omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf      = infile_reader(infile_coh_hf)
        num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)
        
        #------- EEMs Normalization ------
        ## incoherent HF and DFT
        # auc_val_cb_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC before normalization: $auc_val_cb_hf")
        # norm_probs_matrix_incoh_hf = probs_matrix_incoh_hf ./ auc_val_cb_hf
        # # auc_val_ca_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC after normalization: $auc_val_ca_hf")

        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC before normalization: $auc_val_cb_dft")
        norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        # auc_val_ca_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC after normalization: $auc_val_ca_dft")

        # ## coherent HF and DFT
        # auc_val_qb_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf)
        # # println("coherent-HF AUC before normalization: $auc_val_qb_hf")
        # norm_probs_matrix_coh_hf = probs_matrix_coh_hf ./ auc_val_qb_hf
        # # auc_val_qa_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # # println("coherent-HF AUC after normalization: $auc_val_qa_hf")

        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        # println("coherent-DFT AUC before normalization: $auc_val_qb_dft")
        norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        # auc_val_qa_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)
        # println("coherent-DFT AUC after normalization: $auc_val_qa_dft")

        ##------ binning the EEMs ------
        nbins_x  = 200
        nbins_ir = 200
        # binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf          = binning_2d(nbins_x, nbins_ir, omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft       = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)  
        
        ##------ data analysis ------
        n::Int = 10
        # gen_path   = "./binned_EEMs/cdte_$(isys)/peaks_data/"
        # mkpath(gen_path)
        # path_incoh_hf  = gen_path * "cdte$(isys)_incoherent_hf_max_peaks_location.out"
        # path_incoh_dft = gen_path * "cdte$(isys)_incoherent_b3lyp_max_peaks_location.out"
        # path_coh_hf  = gen_path * "cdte$(isys)_coherent_hf_max_peaks_location.out"
        # path_coh_dft = gen_path * "cdte$(isys)_coherent_b3lyp_max_peaks_location.out"

        # writing_max_peaks_location_to_outfile(n, au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf, path_incoh_hf)
        # writing_max_peaks_location_to_outfile(n, au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, path_incoh_dft)
        # writing_max_peaks_location_to_outfile(n, au_to_ev, binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf, path_coh_hf)
        # writing_max_peaks_location_to_outfile(n, au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, path_coh_dft)
        
        ##------ EEMs max peaks extraction ------
        # mat_incoh_hf[idx,:]  = single_max_peak_location(binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf)
        # mat_coh_hf[idx,:]  = single_max_peak_location(binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf)
        # mat_incoh_dft[idx,:] = single_max_peak_location(binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft)
        # mat_coh_dft[idx,:] = single_max_peak_location(binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft)

        ##------ cEEMs Table data extraction ------
        # conditional_mat_incoh_hf  = get_conditional_mat_v2(binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf)
        # conditional_mat_coh_hf  = get_conditional_mat_v2(binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf)
        conditional_mat_incoh_dft = get_conditional_mat_v2(binned_probs_matrix_incoh_dft)
        conditional_mat_coh_dft = get_conditional_mat_v2(binned_probs_matrix_coh_dft)
        
        # mat_incoh_hf[idx,:]  = single_max_peak_location(binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, conditional_mat_incoh_hf)
        # mat_coh_hf[idx,:]  = single_max_peak_location(binned_omega_x_coh_hf, binned_omega_ir_coh_hf, conditional_mat_coh_hf)
        mat_incoh_dft[idx,:] = single_max_peak_location(binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, conditional_mat_incoh_dft)
        mat_coh_dft[idx,:] = single_max_peak_location(binned_omega_x_coh_dft, binned_omega_ir_coh_dft, conditional_mat_coh_dft)

        # ##------ IR marginals max peaks extraction ------
        # marg_ir_incoh_hf, marg_ir_coh_hf   = get_marginal_ir_binned_eem(binned_probs_matrix_incoh_hf, binned_probs_matrix_coh_hf)
        # marg_ir_incoh_dft, marg_ir_coh_dft = get_marginal_ir_binned_eem(binned_probs_matrix_incoh_dft, binned_probs_matrix_coh_dft)
        # mat_incoh_hf[idx,:]  = single_max_peak_location_v2(binned_omega_ir_incoh_hf, marg_ir_incoh_hf)
        # mat_coh_hf[idx,:]  = single_max_peak_location_v2(binned_omega_ir_coh_hf, marg_ir_coh_hf)
        # mat_incoh_dft[idx,:] = single_max_peak_location_v2(binned_omega_ir_incoh_dft, marg_ir_incoh_dft)
        # mat_coh_dft[idx,:] = single_max_peak_location_v2(binned_omega_ir_coh_dft, marg_ir_coh_dft)

    end

    #----- writing EEMs peaks to outfile ------
    gen_dir_path   = "./first_order_corr_res/tabels/"
    # mkpath(gen_dir_path)

    # out_path_hf = gen_dir_path * "EEMs_global_max_hf.out"
    # writing_max_peak_table_v3(mat_incoh_hf, mat_coh_hf, out_path_hf)

    # out_path_dft = gen_dir_path * "EEMs_global_max_b3lyp.out"
    # writing_max_peak_table_v3(mat_incoh_dft, mat_coh_dft, out_path_dft)

    # #----- writing cEEMs peaks to outfile ------
    # gen_dir_path   = "./binned_EEMs/tables_data/"
    # # mkpath(gen_dir_path)

    # out_path_hf = gen_dir_path * "cEEMs_global_max_hf.out"
    # writing_max_peak_table_v4(mat_incoh_hf, mat_coh_hf, out_path_hf)

    out_path_dft = gen_dir_path * "cEEMs_global_max_b3lyp.out"
    writing_max_peak_table_v4(mat_incoh_dft, mat_coh_dft, out_path_dft)

    # #----- writing IR marginals peaks to outfile ------
    # gen_dir_path   = "./binned_EEMs/tables_data/"
    # # mkpath(gen_dir_path)

    # out_path_hf = gen_dir_path * "IR_marginals_global_max_hf.out"
    # writing_max_peak_table_v5(mat_incoh_hf, mat_coh_hf, out_path_hf)

    # out_path_dft = gen_dir_path * "IR_marginals_global_max_b3lyp.out"
    # writing_max_peak_table_v5(mat_incoh_dft, mat_coh_dft, out_path_dft)

    return nothing
end
##================================================================##
function test2A()
##================================================================##
## Binning the EEMs and plotting the unitless EEMs!
    sys_list = [100]
    # gen_dir_path = "./first_order_corr_res/decompos_plt/"
    # mkpath(gen_dir_path)

    nbins_x  = 200
    nbins_ir = 200

    for (idx, isys) in enumerate(sys_list)
        infiles_dir_path = "./elec_hole_corrected/elec_hole_corrected_raw_data/cdte_$(isys)/"

        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        # infile_incoh_hf  = infiles_dir_path * "cdte$(isys)_hf_classical_prob_res_0.0625.out"
        # infile_coh_hf    = infiles_dir_path * "cdte$(isys)_hf_quantum_prob_res_0.0625.out"
        infile_incoh_dft = infiles_dir_path * "cdte$(isys)_b3lyp_classical_prob_res_0.0625.out"
        infile_coh_dft   = infiles_dir_path * "cdte$(isys)_b3lyp_quantum_prob_res_0.0625.out"

        # num_x_incoh_hf, num_ir_incoh_hf, omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf      = infile_reader(infile_incoh_hf)
        # num_x_coh_hf, num_ir_coh_hf, omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf      = infile_reader(infile_coh_hf)
        num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)
        
        #------- EEMs Normalization ------
        # ## incoherent HF and DFT
        # auc_val_cb_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC before normalization: $auc_val_cb_hf")
        # norm_probs_matrix_incoh_hf = probs_matrix_incoh_hf ./ auc_val_cb_hf
        # # auc_val_ca_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC after normalization: $auc_val_ca_hf")

        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC before normalization: $auc_val_cb_dft")
        norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        # auc_val_ca_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC after normalization: $auc_val_ca_dft")

        # ## coherent HF and DFT
        # auc_val_qb_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf)
        # # println("coherent-HF AUC before normalization: $auc_val_qb_hf")
        # norm_probs_matrix_coh_hf = probs_matrix_coh_hf ./ auc_val_qb_hf
        # # auc_val_qa_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # # println("coherent-HF AUC after normalization: $auc_val_qa_hf")

        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        # println("coherent-DFT AUC before normalization: $auc_val_qb_dft")
        norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        # auc_val_qa_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)
        # println("coherent-DFT AUC after normalization: $auc_val_qa_dft")

        ##------ binning the EEMs ------
        # binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)  
        
        #----- surface plotting the binned EEMs ------
        # println("Plotting the binned EEMs for system size: $isys")
        # plt_dir_path   = "./first_order_corr_res/cdte_$(isys)/surf_plts/"
        # # mkpath(plt_dir_path)
        # plt_path_incoh_hf  = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_incoherent_hf.png"
        # plt_path_coh_hf    = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_coherent_hf.png"
        # plt_path_incoh_dft = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_incoherent_b3lyp.png"
        # plt_path_coh_dft   = plt_dir_path * "cdte_$(isys)_3d_binned_EEM_coherent_b3lyp.png"
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf, plt_path_incoh_hf)
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf, plt_path_coh_hf)
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, plt_path_incoh_dft)
        # surface_3d_plot_binned_eem(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, plt_path_coh_dft)
    
        # #----- 2D plotting max peaks of the binned EEMs ------
        println("Plotting decompositions for system size: $isys")
        plt_dir_path = "./first_order_corr_res/cdte_$(isys)/decompos_plt/"
        mkpath(plt_dir_path)
        plt_path_dft_x  = plt_dir_path * "cdte_$(isys)_incoh_decompos_vs_omegax_b3lyp.pdf"
        plt_path_dft_ir = plt_dir_path * "cdte_$(isys)_incoh_decompos_vs_omegair_b3lyp.pdf"
        # plt_path_dft_x  = gen_dir_path * "cdte_$(isys)_coh_decompos_vs_omegax_b3lyp.pdf"
        # plt_path_dft_ir = gen_dir_path * "cdte_$(isys)_coh_decompos_vs_omegair_b3lyp.pdf"
        
        # plt_path_hf_x  = gen_dir_path * "cdte_$(isys)_prob_max_incoh_vs_omegax_hf.pdf"
        # plt_path_hf_ir = gen_dir_path * "cdte_$(isys)_prob_max_incoh_vs_omegair_hf.pdf"
        # plt_signal_decomposition_v1(binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf, binned_probs_matrix_coh_hf, plt_path_hf_ir, plt_path_hf_x)
        # plt_signal_decomposition_v1(binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, binned_probs_matrix_coh_dft, plt_path_dft_ir, plt_path_dft_x)
    
        ###---------- first and zero order together --------

        #----- 2D plotting max peaks of the binned EEMs ------
        # plt_dir_path = "./first_order_corr_res/decompos_plt/cdte_$(isys)/"
        # mkpath(plt_dir_path)
        # plt_path_dft_x  = plt_dir_path * "cdte_$(isys)_incoh_decompos_vs_omegax_b3lyp"
        # plt_path_dft_ir = plt_dir_path * "cdte_$(isys)_incoh_decompos_vs_omegair_b3lyp_extra.pdf"
        # plt_signal_decomposition_v2(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, binned_probs_matrix_coh_dft, plt_path_dft_ir, plt_path_dft_x)
        

        # plt_path_dft_x  = plt_dir_path * "cdte_$(isys)_coh_decompos_vs_omegax_b3lyp"
        # plt_path_dft_ir = plt_dir_path * "cdte_$(isys)_coh_decompos_vs_omegair_b3lyp_extra.pdf"
        # plt_signal_decomposition_v2(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, binned_probs_matrix_coh_dft, plt_path_dft_ir, plt_path_dft_x)
        
    end
    return nothing
end
##================================================================##
function test3A()
##================================================================##
## matrix representation of EEMs and cEEMs!
    sys_list = [100]
    # gen_dir_path = "./matrix_reps/"
    # mkpath(gen_dir_path)
    for (idx, isys) in enumerate(sys_list)
        # gen_dir_path = "./cdte_$(isys)/binned_EEMs/"
        # mkpath(gen_dir_path)

        au_to_ev = 27.211
        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        infile_incoh_hf  = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_hf_incoherent_prob_res_0.0625.out"
        infile_coh_hf  = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_hf_coherent_prob_res_0.0625.out"
        infile_incoh_dft = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_b3lyp_incoherent_prob_res_0.0625.out"
        infile_coh_dft = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_b3lyp_coherent_prob_res_0.0625.out"

        # num_x_incoh_hf, num_ir_incoh_hf, omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf      = infile_reader(infile_incoh_hf)
        # num_x_coh_hf, num_ir_coh_hf, omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf      = infile_reader(infile_coh_hf)
        num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)
        
        #------- EEMs Normalization ------
        ## incoherent HF and DFT
        # auc_val_cb_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC before normalization: $auc_val_cb_hf")
        # norm_probs_matrix_incoh_hf = probs_matrix_incoh_hf ./ auc_val_cb_hf
        # # auc_val_ca_hf = area_under_the_curve_calculator_2D(omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        # # println("incoherent-HF AUC after normalization: $auc_val_ca_hf")

        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC before normalization: $auc_val_cb_dft")
        norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        # auc_val_ca_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # println("incoherent-DFT AUC after normalization: $auc_val_ca_dft")

        ## coherent HF and DFT
        # auc_val_qb_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, probs_matrix_coh_hf)
        # # println("coherent-HF AUC before normalization: $auc_val_qb_hf")
        # norm_probs_matrix_coh_hf = probs_matrix_coh_hf ./ auc_val_qb_hf
        # # auc_val_qa_hf = area_under_the_curve_calculator_2D(omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        # # println("coherent-HF AUC after normalization: $auc_val_qa_hf")

        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        # println("coherent-DFT AUC before normalization: $auc_val_qb_dft")
        norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        # auc_val_qa_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)
        # println("coherent-DFT AUC after normalization: $auc_val_qa_dft")

        ##------ binning the EEMs ------
        nbins_x  = 30
        nbins_ir = 30
        # binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_incoh_hf, omega_ir_incoh_hf, norm_probs_matrix_incoh_hf)
        binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        # binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf    = binning_2d(nbins_x, nbins_ir, omega_x_coh_hf, omega_ir_coh_hf, norm_probs_matrix_coh_hf)
        binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft) 
        
        ##----- heatmap plotting the binned EEMs ------
        println("Plotting the binned EEMs for system size: $isys")
        gen_dir_path = "./matrix_reps/cdte_$(isys)/"
        # mkpath(gen_dir_path)
        # plt_path_incoh_hf  = gen_dir_path * "cdte_$(isys)_2d_incoherent_mat_rep_hf.pdf"
        # plt_path_coh_hf  = gen_dir_path * "cdte_$(isys)_2d_coherent_mat_rep_hf.pdf"
        plt_path_incoh_dft = gen_dir_path * "cdte_$(isys)_2d_incoherent_mat_rep_b3lyp.pdf"
        plt_path_coh_dft = gen_dir_path * "cdte_$(isys)_2d_coherent_mat_repb3lyp.pdf"
        plt_path_i_dft = gen_dir_path * "cdte_$(isys)_2d_interference_mat_rep_b3lyp.pdf"
        # heatmap_plot_matrix_representation_v1(au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf, plt_path_incoh_hf)
        # heatmap_plot_matrix_representation_v1(au_to_ev, binned_omega_x_coh_hf, binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf, plt_path_coh_hf)
        # heatmap_plot_matrix_representation_v1(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft, plt_path_incoh_dft)
        # heatmap_plot_matrix_representation_v1(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft, plt_path_coh_dft)
        
        # interf_mat_dft = binned_probs_matrix_coh_dft .- binned_probs_matrix_incoh_dft
        # heatmap_plot_matrix_representation_v1(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, interf_mat_dft, plt_path_i_dft)
    
        
        ##-------- cEEMs ------------
        # conditional_mat_incoh_hf  = get_conditional_mat_v2(binned_omega_ir_incoh_hf, binned_probs_matrix_incoh_hf)
        # conditional_mat_coh_hf  = get_conditional_mat_v2(binned_omega_ir_coh_hf, binned_probs_matrix_coh_hf)
        conditional_mat_incoh_dft = get_conditional_mat_v2(binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft)
        conditional_mat_coh_dft = get_conditional_mat_v2(binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft)
        # plt_path_incoh_hf  = gen_dir_path * "cdte_$(isys)_2d_incoherent_cmat_rep_hf.pdf"
        # plt_path_coh_hf  = gen_dir_path * "cdte_$(isys)_2d_coherent_cmat_rep_hf.pdf"
        plt_path_incoh_dft = gen_dir_path * "cdte_$(isys)_2d_incoherent_cmat_rep_b3lyp.pdf"
        plt_path_coh_dft = gen_dir_path * "cdte_$(isys)_2d_coherent_cmat_repb3lyp.pdf"
        plt_path_i_dft = gen_dir_path * "cdte_$(isys)_2d_interference_cmat_rep_b3lyp.pdf"
        
        # val_q, idx_max_q = findmax(conditional_mat_coh_dft)
        # println()
        # println("Maximum value in conditional_mat_coh_dft: $val_q")
        # println("Corresponding omega_x for max in conditional_mat_coh_dft: $(binned_omega_x_coh_dft[idx_max_q[1]] * au_to_ev) eV")
        # println("Corresponding omega_ir for max in conditional_mat_coh_dft: $(binned_omega_ir_coh_dft[idx_max_q[2]] * au_to_ev) eV")
        
        # val_c, idx_max_c = findmax(conditional_mat_incoh_dft) 
        # println()  
        # println("Maximum value in conditional_mat_incoh_dft: $val_c")
        # println("Corresponding omega_x for max in conditional_mat_incoh_dft: $(binned_omega_x_incoh_dft[idx_max_c[1]] * au_to_ev) eV")
        # println("Corresponding omega_ir for max in conditional_mat_incoh_dft: $(binned_omega_ir_incoh_dft[idx_max_c[2]] * au_to_ev) eV")
        # println()  
        
        # prob_range = (0.0, maximum(vcat(vec(conditional_mat_incoh_dft), vec(conditional_mat_coh_dft))))
        # prob_range = (0.0, quantile(vcat(vec(conditional_mat_incoh_dft), vec(conditional_mat_coh_dft)), 0.99))
        # heatmap_plot_matrix_representation_v2(au_to_ev, binned_omega_x_incoh_hf, binned_omega_ir_incoh_hf, conditional_mat_incoh_hf, plt_path_incoh_hf)
        # heatmap_plot_matrix_representation_v2(au_to_ev, binned_omega_x_coh_hf, binned_omega_ir_coh_hf, conditional_mat_coh_hf, plt_path_coh_hf)
        # heatmap_plot_matrix_representation_v2(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, conditional_mat_incoh_dft, plt_path_incoh_dft, prob_range)
        # heatmap_plot_matrix_representation_v2(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, conditional_mat_coh_dft, plt_path_coh_dft, prob_range)
        # heatmap_plot_matrix_representation_v2_gamma(au_to_ev, binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, conditional_mat_incoh_dft, plt_path_incoh_dft, prob_range)
        # heatmap_plot_matrix_representation_v2_gamma(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, conditional_mat_coh_dft, plt_path_coh_dft, prob_range)

        conditional_interf_mat_dft = conditional_mat_coh_dft .- conditional_mat_incoh_dft
        # # interf_absmax = maximum(abs.(conditional_interf_mat_dft))
        # # interf_absmax = quantile(abs.(vec(conditional_interf_mat_dft)), 0.99)
        # # interf_range = (-interf_absmax, interf_absmax)
        # interf_range = (minimum(conditional_interf_mat_dft), maximum(conditional_mat_coh_dft))
        # heatmap_plot_matrix_representation_v3(au_to_ev, binned_omega_x_coh_dft, binned_omega_ir_coh_dft, conditional_interf_mat_dft, plt_path_i_dft, interf_range)

        @show minimum(conditional_interf_mat_dft)
        @show maximum(conditional_interf_mat_dft)
        # val_i, idx_max_i = findmax(conditional_interf_mat_dft) 
        # println()  
        # println("Maximum value in conditional_interf_mat_dft: $val_i")
        # println("Corresponding omega_x for max in conditional_interf_mat_dft: $(binned_omega_x_coh_dft[idx_max_i[1]] * au_to_ev) eV")
        # println("Corresponding omega_ir for max in conditional_interf_mat_dft: $(binned_omega_ir_coh_dft[idx_max_i[2]] * au_to_ev) eV")

    end
    return nothing
end
##================================================================##
function test4A()
##================================================================##
## Binning the EEMs and plotting the unitless EEMs!
    sys_list = [20]
    # gen_dir_path = "./first_order_corr_res/decompos_plt/"
    # mkpath(gen_dir_path)

    nbins_x  = 200
    nbins_ir = 200

    for (idx, isys) in enumerate(sys_list)
        infiles_dir_path = "./elec_hole_corrected/elec_hole_corrected_raw_data/cdte_$(isys)/"

        ##----- Reading the HF and DFT raw data for incoherent and coherent -------
        # infile_incoh_hf  = infiles_dir_path * "cdte$(isys)_hf_classical_prob_res_0.0625.out"
        # infile_coh_hf    = infiles_dir_path * "cdte$(isys)_hf_quantum_prob_res_0.0625.out"
        infile_incoh_dft = infiles_dir_path * "cdte$(isys)_b3lyp_classical_prob_res_0.0625.out"
        infile_coh_dft   = infiles_dir_path * "cdte$(isys)_b3lyp_quantum_prob_res_0.0625.out"

        infile_incoh_zo_dft = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_b3lyp_classical_prob_res_0.0625.out"
        infile_coh_zo_dft = "./transition_prob_eems/cdte_$(isys)/raw_data/cdte$(isys)_b3lyp_quantum_prob_res_0.0625.out"


        num_x_incoh_dft, num_ir_incoh_dft, omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft = infile_reader(infile_incoh_dft)
        num_x_coh_dft, num_ir_coh_dft, omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft = infile_reader(infile_coh_dft)
        num_x_incoh_zo_dft, num_ir_incoh_zo_dft, omega_x_incoh_zo_dft, omega_ir_incoh_zo_dft, probs_matrix_incoh_zo_dft = infile_reader(infile_incoh_zo_dft)
        num_x_coh_zo_dft, num_ir_coh_zo_dft, omega_x_coh_zo_dft, omega_ir_coh_zo_dft, probs_matrix_coh_zo_dft = infile_reader(infile_coh_zo_dft)
        #------- EEMs Normalization ------
        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_dft, omega_ir_incoh_dft, probs_matrix_incoh_dft)
        norm_probs_matrix_incoh_dft = probs_matrix_incoh_dft ./ auc_val_cb_dft
        
        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_dft, omega_ir_coh_dft, probs_matrix_coh_dft)
        norm_probs_matrix_coh_dft = probs_matrix_coh_dft ./ auc_val_qb_dft
        
        auc_val_cb_dft = area_under_the_curve_calculator_2D(omega_x_incoh_zo_dft, omega_ir_incoh_zo_dft, probs_matrix_incoh_zo_dft)
        norm_probs_matrix_incoh_zo_dft = probs_matrix_incoh_zo_dft ./ auc_val_cb_dft
        
        auc_val_qb_dft = area_under_the_curve_calculator_2D(omega_x_coh_zo_dft, omega_ir_coh_zo_dft, probs_matrix_coh_zo_dft)
        norm_probs_matrix_coh_zo_dft = probs_matrix_coh_zo_dft ./ auc_val_qb_dft
        
        ##------ binning the EEMs ------
        binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_dft, omega_ir_incoh_dft, norm_probs_matrix_incoh_dft)
        binned_omega_x_coh_dft, binned_omega_ir_coh_dft, binned_probs_matrix_coh_dft = binning_2d(nbins_x, nbins_ir, omega_x_coh_dft, omega_ir_coh_dft, norm_probs_matrix_coh_dft)  
        binned_omega_x_incoh_zo_dft, binned_omega_ir_incoh_zo_dft, binned_probs_matrix_incoh_zo_dft = binning_2d(nbins_x, nbins_ir, omega_x_incoh_zo_dft, omega_ir_incoh_zo_dft, norm_probs_matrix_incoh_zo_dft)
        binned_omega_x_coh_zo_dft, binned_omega_ir_coh_zo_dft, binned_probs_matrix_coh_zo_dft = binning_2d(nbins_x, nbins_ir, omega_x_coh_zo_dft, omega_ir_coh_zo_dft, norm_probs_matrix_coh_zo_dft)  
      
        ###---------- first and zero order together --------
        println("Plotting decompositions for system size: $isys")
        plt_dir_path = "./first_order_corr_res/cdte_$(isys)/all/"
        mkpath(plt_dir_path)
        plt_path_dft_x = plt_dir_path * "cdte_$(isys)_incoh_decompos_vs_omegax_b3lyp.pdf"
       
        plt_signal_decomposition_v3(binned_omega_x_incoh_dft, binned_omega_ir_incoh_dft, binned_probs_matrix_incoh_zo_dft, binned_probs_matrix_coh_zo_dft, binned_probs_matrix_incoh_dft, binned_probs_matrix_coh_dft, plt_path_dft_x)
    
        
        
    end
    return nothing
end
##================================================================##
## MAIN CODE ##
##================================================================##
println("start_time= ",now())
println()
    test4A()
println()
println("end_time= ",now())
