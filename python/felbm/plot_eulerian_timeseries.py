from analyze_eulerian_timeseries import *
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Plot pos")
    parser.add_argument("folder", type=str, help="Folder")
    parser.add_argument("--save", action="store_true", help="Save figures")
    parser.add_argument("--show", action="store_true", help="Show figures")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()

    felbm_folder = args.folder
    analysisfolder = os.path.join(felbm_folder, "Analysis")

    if rank==0:
        if True:
            tdata = np.loadtxt(os.path.join(analysisfolder, "tdata.dat"))
            t_ = tdata[0, :]
            uxt = tdata[1, :]
            uyt = tdata[2, :]
            S1t = tdata[3, :]
            S2t = tdata[4, :]

            with h5py.File(os.path.join(analysisfolder, "pore_data.h5"), "r") as h5f:
                un = np.array(h5f["un"])
                n = np.array(h5f["n"])
                t = np.array(h5f["t"])
                lp = np.array(h5f["lp"])
                d = h5f.attrs.get("d")
            #print(un.shape)
            print(lp)
            #for i in range(len(un)// 10000):
            uy_mean = uyt.mean()
            print("uy_mean = {}".format(uy_mean))
            print("d       = {}".format(d))
            t_adv = d / uy_mean
            print("t_adv   = {}".format(t_adv))

            # exit()

            unmean = un.mean(axis=1)
            dun = un - np.outer(unmean, np.ones_like(t))
            import time
            from scipy.interpolate import InterpolatedUnivariateSpline

            do_plot = False
            tol = 1e-14

            ptime0 = time.time()
            data = [None for _ in range(len(un))]
            valid = [False for _ in range(len(un))]
            
            for i in range(len(un)):
                if any(abs(dun[i, :-1] - dun[i, 1:]) < tol):
                    continue
                valid[i] = True
                #if True:
                dup = dun[i, :] > 0
                change = dup[1:] ^ dup[:-1]
                t_x = (t[1:] * dun[i, :-1] - t[:-1] * dun[i, 1:]) / (dun[i, :-1] - dun[i, 1:])
                t_x = t_x[change]
                tau = t_x[1:] - t_x[:-1]

                fintp = InterpolatedUnivariateSpline(t, dun[i, :], k=1)
                A = np.zeros_like(tau)
                for k in range(len(tau)):
                    A[k] = fintp.integral(t_x[k], t_x[k+1])
                    #print(A)
                data[i] = np.vstack((A, tau)).T

                if i < 5:
                    fig, ax = plt.subplots(1, 1)
                    ax.plot(t, dun[i, :])
                    ax.plot(t, 0*t)
                    for k in range(len(tau)):
                        un_x = A[k]/tau[k]
                        ax.plot([t_x[k], t_x[k+1]], [un_x, un_x])
            
                    if args.save:
                        plt.savefig(os.path.join(analysisfolder, "example_{}_network.png".format(i)))

                    if args.show:
                        plt.show()

            tau_ = []
            A_ = []

            tau_mean = np.zeros(len(data))
            A_mean = np.zeros_like(tau_mean)
            dun_mean = np.zeros_like(tau_mean)
            for i, Atau in enumerate(data):
                if valid[i] and Atau.shape[0] > 0:
                    #print(Atau.shape)
                    tau_.extend(Atau[:, 1].flatten())
                    A_.extend(Atau[:, 0].flatten())
                    tau_mean[i] = Atau[:, 1].mean()
                    A_mean[i] = Atau[:, 0].mean()
                    dun_mean[i] = (Atau[:, 0]/Atau[:, 1]).mean()
            theta = np.arcsin((n[valid, 1]))

            A_ = np.array(A_)
            tau_ = np.array(tau_)

            tau_avg = tau_.mean()
            A_avg = abs(A_).mean()
            dun_avg = (abs(A_)/tau_).mean()

            print("tau_avg  = {}".format(tau_avg))
            print("A_avg    = {}".format(A_avg))
            print("dun_avg  = {}".format(dun_avg))

            if True:
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
                ax1.plot(theta, unmean[valid], ',')
                ax1.set_xlabel("$\\theta$")
                ax1.set_ylabel("$\\bar{u}_n$")
                ax2.plot(theta, dun_mean[valid], ',')
                ax2.set_xlabel("$\\theta$")
                ax2.set_ylabel("$\\Delta u_n$")
                ax3.plot(theta, tau_mean[valid], ',')
                ax3.set_xlabel("$\\theta$")
                ax3.set_ylabel("$\\bar{\\tau}$")
                ax4.plot(theta, A_mean[valid], ',')
                ax4.set_xlabel("$\\theta$")
                ax4.set_ylabel("$\\bar{A}$")

                if args.save:
                    plt.savefig(os.path.join(analysisfolder, "angle_dependency_network.png"))

                if args.show:
                    plt.show()

                fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
                ax1.hist(tau_, bins=100, density=True)
                ax1.set_yscale("log")
                ax1.set_xlabel("$\\tau$")
                ax1.set_ylabel("$P(\\tau)$")
                x = np.linspace(0., tau_.max(), 1000)
                ax1.plot(x, 1./tau_avg * np.exp(-x/tau_avg))

                ax2.hist(abs(A_), bins=100, density=True)
                ax2.set_yscale("log")
                ax2.set_xlabel("$A$")
                ax2.set_ylabel("$P(A)$")
                x = np.linspace(0., abs(A_).max(), 1000)
                ax2.plot(x, 1./A_avg * np.exp(-x/A_avg))
                
                ax3.hist(abs(A_)/tau_, bins=100, density=True)
                ax3.set_yscale("log")
                ax3.set_xlabel("$|\Delta u_n|$")
                ax3.set_ylabel("$P(|\Delta u_n|)$")
                x = np.linspace(0., (A_/tau_).max(), 1000)
                ax3.plot(x, 1./dun_avg * np.exp(-x/dun_avg))

                if args.save:
                    plt.savefig(os.path.join(analysisfolder, "Pdfs_network.png"))

                if args.show:
                    plt.show()
                
            nbins = 300
            dtau = t_adv/100
            As = [[] for _ in range(nbins)]
            uns = [[] for _ in range(nbins)]
            for Ai, taui in zip(A_, tau_):
                i = int(np.round(taui/dtau))
                if i < len(As):
                    As[i].append(abs(Ai))
                    uns[i].append(abs(Ai/taui))
            As_mean = []
            uns_mean = []
            for i in range(nbins):
                if len(As[i]) == 0:
                    break
                As_mean.append(np.mean(As[i]))
                uns_mean.append(np.mean(uns[i]))
            As_mean = np.array(As_mean)
            uns_mean = np.array(uns_mean)
            taus = dtau/2 + dtau * np.arange(len(As_mean))

            def fitfunc(xx, u0, tau0):
                return u0 * (1.0 - np.exp(-xx/tau0)) 
            dun_ = (abs(A_)/tau_)

            from scipy.optimize import curve_fit
            popt, pcov = curve_fit(fitfunc, taus, uns_mean, p0=[1., t_adv])

            print("fitted parameters:")
            print("u0   = {}".format(popt[0]))
            print("tau0 = {}".format(popt[1]))

            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.plot(tau_, (abs(A_)), ',')
            ax1.plot(taus, As_mean)
            ax1.plot(taus, taus * fitfunc(taus, *popt))
            ax1.set_xlabel("$\\tau$")
            ax1.set_ylabel("$A$")

            ax2.plot((tau_), dun_, ',')
            ax2.plot(taus, uns_mean)
            ax2.plot(taus, fitfunc(taus, *popt))
            ax2.set_xlabel("$\\tau$")
            ax2.set_ylabel("$\\Delta u_n$")

            if args.save:
                plt.savefig(os.path.join(analysisfolder, "scatterfit_network.png"))

            if args.show:
                plt.show()


    if False:
        if True:
            print("Time --> Space:")

            P2f = np.loadtxt(os.path.join(analysisfolder, "P2f.dat"))
            fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            df = P2f[1, 0]-P2f[0, 0]
            Pf = np.sqrt(P2f[:, 1])
            Pf /= Pf.sum() * df
            Pfy = np.sqrt(P2f[:, 2])
            Pfy /= Pfy.sum() * df
            for i, ax in enumerate(axs):
                ax.plot(P2f[:, 0], Pf, label="")
                ax.plot(P2f[:, 0], Pfy, label="$u_y$")
                ax.set_xlabel("f")
                ax.set_ylabel("P(f)")
                if i in [1, 2]:
                    ax.set_yscale("log")
                if i in [2]:
                    ax.set_xscale("log")
                ax.legend()

            if args.save:
                plt.savefig(os.path.join(analysisfolder, "Pf.png"))

            if args.show:
                plt.show()

            ux = np.loadtxt(os.path.join(analysisfolder, "uxnorm_avg.dat"))
            uy = np.loadtxt(os.path.join(analysisfolder, "uynorm_avg.dat"))
            freq = np.loadtxt(os.path.join(analysisfolder, "freq_avg.dat"))
            freqy = np.loadtxt(os.path.join(analysisfolder, "freqy_avg.dat"))
            rho = np.loadtxt(os.path.join(analysisfolder, "rho_avg.dat"))
            du = np.loadtxt(os.path.join(analysisfolder, "du.dat"))
            is_fluid_xy = np.array(np.loadtxt(os.path.join(analysisfolder, "is_fluid_xy.dat")), dtype=bool)

            tol = 1e-9

            Ax = np.zeros_like(ux)
            Ax[ux > tol] = ux[ux > tol]/freq[ux > tol]

            Ay = np.zeros_like(uy)
            Ay[uy > tol] = uy[uy > tol]/freqy[uy > tol]

            freq_mean = freq[is_fluid_xy].mean()
            freqy_mean = freqy[is_fluid_xy].mean()
            ux_mean = ux[is_fluid_xy].mean()
            Ax_mean = Ax[is_fluid_xy][freq[is_fluid_xy] > tol].mean()
            uy_mean = uy[is_fluid_xy].mean()
            Ay_mean = Ay[is_fluid_xy][freqy[is_fluid_xy] > tol].mean()
            du_mean = du[is_fluid_xy].mean()
            print("freq_mean =", freq_mean)
            print("ux_mean =", ux_mean)
            print("Ax_mean =", Ax_mean)
            print("freqy_mean =", freqy_mean)
            print("uy_mean =", uy_mean)
            print("Ay_mean =", Ay_mean)
            print("du_mean =", du_mean)

            fig, ((ax0, ax1, ax2, ax3), (ax4, ax5, ax6, ax7)) = plt.subplots(2, 4, figsize=(20, 10))

            ax0.hist(freq[is_fluid_xy], bins=256, density=True)
            ax0.axvline(freq_mean, c="k")
            #ax0.set_yscale("log")
            ax0.set_xlabel("$\\bar{f}$")
            ax0.set_ylabel("$P(\\bar{f})$")

            ax1.hist(freqy[is_fluid_xy], bins=256, density=True)
            ax1.axvline(freqy_mean, c="k")
            #ax1.set_yscale("log")
            ax1.set_xlabel("$\\bar{f}$ ($u_y$)")
            ax1.set_ylabel("$P(\\bar{f})$ ($u_y$)")

            ax2.hist(ux[is_fluid_xy], bins=256, density=True)
            ax2.axvline(ux_mean, c="k")
            ax2.set_yscale("log")
            ax2.set_xlabel("$\\bar{u}_x$")
            ax2.set_ylabel("$P(\\bar{u}_x)$")
            x_ = np.linspace(0., 10*ux_mean, 100)
            ax2.plot(x_, 1./ux_mean*np.exp(-x_/ux_mean))

            ax3.hist(Ax[is_fluid_xy][freq[is_fluid_xy] > tol], bins=256, density=True)
            ax3.axvline(Ax_mean, c="k")
            x_ = np.linspace(0., 10*Ax_mean, 100)
            ax3.plot(x_, 1./Ax_mean*np.exp(-x_/Ax_mean))
            ax3.set_xlabel("$A_x$")
            ax3.set_ylabel("$P(A_x)$")
            ax3.set_yscale("log")

            ax4.hist(uy[is_fluid_xy], bins=256, density=True)
            ax4.axvline(uy_mean, c="k")
            ax4.set_yscale("log")
            ax4.set_xlabel("$\\bar{u}_y$")
            ax4.set_ylabel("$P(\\bar{u}_y)$")
            x_ = np.linspace(0., 10*uy_mean, 100)
            ax4.plot(x_, 1./uy_mean*np.exp(-x_/uy_mean))

            ax5.hist(Ay[is_fluid_xy][freqy[is_fluid_xy] > tol], bins=256, density=True)
            ax5.axvline(Ay_mean, c="k")
            x_ = np.linspace(0., 10*Ay_mean, 100)
            ax5.plot(x_, 1./Ay_mean*np.exp(-x_/Ay_mean))
            ax5.set_xlabel("$A_y$")
            ax5.set_ylabel("$P(A_y)$")
            ax5.set_yscale("log")

            ax6.hist(du[is_fluid_xy], bins=256, density=True)
            ax6.axvline(du_mean, c="k")
            #x_ = np.linspace(0., 10*du_mean, 100)
            #ax4.plot(x_, 1./du_mean * np.exp(-x_/du_mean))
            ax6.set_xlabel("$u'$")
            ax6.set_ylabel("$P(u')$")
            #ax4.set_yscale("log")

            if args.save:
                plt.savefig(os.path.join(analysisfolder, "pdfs.png"))

            if args.show:
                plt.show()

            #plt.imshow(ux.reshape((ny, nx)))
            fig1, ax1 = plt.subplots(1, 1, figsize=(20,6))

            im1 = ax1.imshow(freq)
            fig1.colorbar(im1, ax=ax1)
            ax1.set_title("freq")

            if args.save:
                plt.savefig(os.path.join(analysisfolder, "spatial_freq.png"))

            if args.show:
                plt.show()

            fig2, ((ax2, ax3), (ax4, ax5)) = plt.subplots(2, 2, figsize=(20,6))
            
            im2 = ax2.imshow(ux)
            ax2.set_title("ux")
            fig2.colorbar(im2, ax=ax2)
            
            im3 = ax3.imshow(Ax)
            ax3.set_title("Ax")
            fig2.colorbar(im3, ax=ax3)
            
            im4 = ax4.imshow(rho)
            ax4.set_title("rho")
            fig2.colorbar(im4, ax=ax4)
            
            im5 = ax5.imshow(du)
            ax5.set_title("du")
            fig2.colorbar(im5, ax=ax5)
            
            if args.save:
                plt.savefig(os.path.join(analysisfolder, "spatial.png"))

            if args.show:
                plt.show()

        #except:
        #    pass

        if True:
            print("Space --> Time:")

            ux_mean = uxt.mean()
            uy_mean = uyt.mean()
            S1_mean = S1t.mean()
            S2_mean = S2t.mean()
            print("ux_avg =", ux_mean)
            print("uy_avg =", uy_mean)
            print("S1_avg =", S1_mean)
            print("S2_avg =", S2_mean)

            fig, (ax1, ax2) = plt.subplots(2, 1)
            ax1.set_ylabel("Mean velocity")
            ax1.set_xlabel("Time t")
            ax1.plot(t_, uxt, label="ux")
            ax1.plot(t_, ux_mean * np.ones_like(t_), label="ux_mean")
            ax1.plot(t_, uyt, label="uy")
            ax1.plot(t_, uy_mean * np.ones_like(t_), label="uy_mean")
            ax1.legend()

            ax2.set_ylabel("Mean cluster size")
            ax2.set_xlabel("Time t")
            ax2.plot(t_, S1t, label="phase1")
            ax2.plot(t_, S1_mean * np.ones_like(t_), label="phase1_mean")
            ax2.plot(t_, S2t, label="phase2")
            ax2.plot(t_, S2_mean * np.ones_like(t_), label="phase2_mean")
            ax2.legend()

            if args.save:
                plt.savefig(os.path.join(analysisfolder, "time_signals.png"))

            if args.show:
                plt.show()

        # except:
        #     pass