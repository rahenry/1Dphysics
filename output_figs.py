#!/usr/bin/env python3


from output import *

def fig_XYW_1(S):
    matplotlib.use("pgf")
    def graph(axs, X, Y, eta, kappa):
        ax = axs[Y][X]
        s = S.find({"eta": eta, "kappa": kappa})[0]
        L = s["L"]
        xdata_sim = [L * x[0].real for x in s["eigs"]]
        ydata_sim = [L * x[0].imag for x in s["eigs"]]
        xdata_exact = [L * x.real for x in s.spectrum]
        ydata_exact = [L * x.imag for x in s.spectrum]
        plot_args = dict(linestyle="", linewidth=10, markeredgewidth=1)
        ax.plot(
            xdata_exact, ydata_exact, marker="x", markersize=8, color="red", **plot_args
        )
        ax.plot(
            xdata_sim, ydata_sim, marker="o", markersize=2, color="blue", **plot_args
        )
        # ax.set_xlabel(r"$\mathbb{R}(E)$")
        xlab = rf"$\Re(E)$"
        info = rf"$\eta = {eta:1.1f}$, $\kappa={kappa:1.1f}$"
        ax.set_title(info)
        if Y > 0:
            ax.set_xlabel(xlab)
        if X == 0:
            ax.set_ylabel("$\Im(E)$")
        ax.axis("equal")
        caption = rf"$\eta = {eta:1.1f}$, $\kappa={kappa:1.1f}$"
        return caption

    fig, axs = plt.subplots(2, 2)
    c = S.config_systems
    caption = ""
    a = graph(axs, 0, 0, c["eta"][1], c["kappa"][1])
    caption += rf"{a} (top left), "
    a = graph(axs, 0, 1, c["eta"][2], c["kappa"][0])
    caption += rf"{a} (top right), "
    a = graph(axs, 1, 0, c["eta"][0], c["kappa"][2])
    caption += rf"{a} (bottom left), "
    a = graph(axs, 1, 1, c["eta"][1], c["kappa"][0])
    caption += rf"{a} (bottom right)"
    plt.tight_layout()
    # plt.show()
    f = os.path.join(S.graph_dir, f"fig_XYW.pdf")
    plt.savefig(f)
    plt.savefig(f.replace(".pdf", ".pgf"))

    caption = rf"Spectrum for the $XYW$ model with various parameters. Blue dots show exact diagonalisation values, and red crosses show the values predicted from the polynomial expressions."
    caption = r"\caption{" + caption + r"}"
    open(f.replace(".pdf", "_caption.tex"), "w+").write(caption)
    # tikzplotlib.save(f.replace(".pdf", ".tikz"))

    plt.close()


def fig_XYW_periodic(S):
    matplotlib.use("pgf")
    def graph(ax, s):
        kappa = s['kappa']
        eta = s['eta']
        bc = s['bc']
        # s = S.find({"eta": eta, "kappa": kappa, "bc":bc})[0]
        L = s["L"]
        xdata_sim = [L * x[0].real for x in s["eigs"]]
        ydata_sim = [L * x[0].imag for x in s["eigs"]]
        xdata_exact = [L * x.real for x in s.spectrum]
        ydata_exact = [L * x.imag for x in s.spectrum]
        plot_args = dict(linestyle="", linewidth=10, markeredgewidth=1)
        ax.plot(
            xdata_exact, ydata_exact, marker="x", markersize=8, color="red", **plot_args
        )
        ax.plot(
            xdata_sim, ydata_sim, marker="o", markersize=2, color="blue", **plot_args
        )
        # ax.set_xlabel(r"$\mathbb{R}(E)$")
        xlab = rf"$\Re(E)$"
        info = rf"$\eta = {eta:1.1f}$, $\kappa={kappa:1.1f}$, bc$={bc}$"
        # ax.set_title(info)
        # if Y > 0:
        #     ax.set_xlabel(xlab)
        # if X == 0:
        #     ax.set_ylabel("$\Im(E)$")
        ax.set_title(f'L={s["L"]}')
        ax.axis("equal")
        caption = rf"$\eta = {eta:1.1f}$, $\kappa={kappa:1.1f}$, bc$={bc}$"
        return caption

    fig, axs = plt.subplots(2, 2)
    c = S.config_systems
    caption = ""
    a = graph(axs[0][0], S.systems[0])
    a = graph(axs[1][0], S.systems[1])
    a = graph(axs[0][1], S.systems[2])
    a = graph(axs[1][1], S.systems[3])
    # a = graph(axs, 0, 0, c["eta"][1], c["kappa"][0], 0)
    # caption += rf"{a} (top left), "
    # a = graph(axs, 0, 1, c["eta"][0], c["kappa"][1], 0)
    # caption += rf"{a} (top right), "
    # a = graph(axs, 1, 0, c["eta"][1], c["kappa"][0], 1)
    # caption += rf"{a} (bottom left), "
    # a = graph(axs, 1, 1, c["eta"][0], c["kappa"][1], 1)
    # caption += rf"{a} (bottom right)"
    # plt.show()
    f = os.path.join(S.graph_dir, f"fig_XYW.pdf")
    print(f)
    xlab = rf"$\Re(E)$"
    ylab = rf"$\Im(E)$"
    axs[1][0].set_ylabel(ylab)
    axs[0][0].set_ylabel(ylab)
    axs[1][1].set_xlabel(xlab)
    axs[1][0].set_xlabel(xlab)
    plt.tight_layout()
    plt.savefig(f)
    plt.savefig(f.replace(".pdf", ".pgf"))

    caption = rf"Spectrum for the $XYW$ model with various parameters. Blue dots show exact diagonalisation values, and red crosses show the values predicted from the polynomial expressions."
    caption = r"\caption{" + caption + r"}"
    open(f.replace(".pdf", "_caption.tex"), "w+").write(caption)
    # tikzplotlib.save(f.replace(".pdf", ".tikz"))

    plt.close()
