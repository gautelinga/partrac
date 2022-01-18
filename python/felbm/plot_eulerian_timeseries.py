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

            tdata = np.loadtxt(os.path.join(analysisfolder, "tdata.dat"))
            t_ = tdata[0, :]
            uxt = tdata[1, :]
            uyt = tdata[2, :]
            S1t = tdata[3, :]
            S2t = tdata[4, :]

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