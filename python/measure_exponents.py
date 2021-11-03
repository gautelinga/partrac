import argparse
import os

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(
        description="Measure exponents from elongation data")
    parser.add_argument("folder", type=str, help="Folder")
    parser.add_argument("-id0",
                        type=int,
                        default=1,
                        help="First id to include")
    parser.add_argument("--show", action="store_true", help="Show")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    filename = os.path.join(args.folder, "Analysis", "elongdata.dat")
    d = np.loadtxt(filename)

    t = d[args.id0:, 0]
    elong_mean = d[args.id0:, 1]
    elong_var = d[args.id0:, 2]

    logelong_mean = d[args.id0:, 4]
    logelong_var = d[args.id0:, 5]

    # print(d)
    idt = (len(t) - args.id0) // 2

    t_fit = t[idt:]
    p_lm = np.polyfit(t_fit, logelong_mean[idt:], 1)
    p_lv = np.polyfit(t_fit, logelong_var[idt:], 1)
    p_em = np.polyfit(t_fit, np.log(elong_mean[idt:]), 1)
    p_ev = np.polyfit(t[idt:], np.log(elong_var[idt:]), 1)

    lam = p_lm[0]
    lam_v = p_lv[0]
    lam_rho = p_em[0]
    lam_rho_v = p_ev[0]

    print("Direct measurements:")
    print("lambda      =", lam)
    print("lambda'     =", lam_v)
    print("lambda_rho  =", lam_rho)
    print("lambda'_rho =", lam_rho_v)

    print("Indirect computations (assuming lognormal rho):")
    lam_ind = 2 * lam_rho - 0.5 * lam_rho_v
    lam_v_ind = -2 * lam_rho + lam_rho_v
    print("lambda  =", lam_ind)
    print("lambda' =", lam_v_ind)

    if args.show:
        fig, ax = plt.subplots(1, 1)

        plt.plot(t, np.log(elong_mean), linestyle="-", label='log(<rho>)')
        plt.plot(t, np.log(elong_var), linestyle="-", label='log(Var(rho))')
        plt.plot(t, logelong_mean, linestyle="-", label='<log(rho)>')
        plt.plot(t, logelong_var, linestyle="-", label='Var(log(rho))')
        plt.plot(t_fit, p_lm[0] * t_fit + p_lm[1], linestyle=":", color="k")
        plt.plot(t_fit, p_lv[0] * t_fit + p_lv[1], linestyle=":", color="k")
        plt.plot(t_fit, p_em[0] * t_fit + p_em[1], linestyle=":", color="k")
        plt.plot(t_fit, p_ev[0] * t_fit + p_ev[1], linestyle=":", color="k")

        plt.legend()
        plt.show()
