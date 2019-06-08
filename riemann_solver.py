#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import json
import math
from scipy import optimize

class Wave:
    def __init__(self, velocity, density, pressure, gamma):
        self.vel = velocity
        self.dens = density
        self.pres = pressure
        self.sonic = math.sqrt(gamma * self.pres / self.dens)

    def __str__(self):
        return "velocity = {0:.3e}, density = {1:.3e}, pressure = {2:.3e}".format(self.vel, self.dens, self.pres)

class RiemannSolver:
    def __init__(self, left, right, gamma):
        self.left = left
        self.right = right
        self.gamma = gamma

    def func(self, p, w):
        if p > w.pres:
            a = 2.0 / ((self.gamma + 1.0) * w.dens)
            b = (self.gamma - 1.0) / (self.gamma + 1.0) * w.pres
            return (p - w.pres) * math.sqrt(a / (p + b))
        else:
            return 2.0 * w.sonic / (self.gamma - 1.0) * (math.pow(p / w.pres, (self.gamma - 1.0) / (2.0 * self.gamma)) - 1.0)

    def solve(self):
        # Two-Rarefaction approximation
        p_init = math.pow((self.left.sonic + self.right.sonic - 0.5 * (self.gamma - 1.0) * (self.right.vel - self.left.vel)) /
                (self.left.sonic / math.pow(self.left.pres, (self.gamma - 1.0) / (2.0 * self.gamma)) +
                self.right.sonic / math.pow(self.right.pres, (self.gamma - 1.0) / (2.0 * self.gamma))), 2.0 * self.gamma / (self.gamma - 1.0))
        f = (lambda p: self.func(p, self.left) + self.func(p, self.right) + self.right.vel - self.left.vel)
        self.pstar = optimize.newton(f, p_init)
        self.vstar = 0.5 * (self.left.vel + self.right.vel) + 0.5 * (self.func(self.pstar, self.right) - self.func(self.pstar, self.left))
        print("f = {0}".format(f(self.pstar)))
        print("pstar = {0}".format(self.pstar))
        print("vstar = {0}".format(self.vstar))
        self.left_shock  = self.pstar > self.left.pres
        self.right_shock = self.pstar > self.right.pres
        if self.left_shock:
            self.densl = self.left.dens * ((self.gamma - 1.0) * self.left.pres + (self.gamma + 1.0)
                    * self.pstar) / ((self.gamma + 1.0) * self.left.pres + (self.gamma - 1.0) * self.pstar)
            self.sl = self.left.vel - self.left.sonic * math.sqrt((self.gamma + 1.0) / (2.0 * self.gamma)
                    * self.pstar / self.left.pres + (self.gamma - 1.0) / (2.0 * self.gamma))
            self.wlstar = Wave(self.vstar, self.densl, self.pstar, self.gamma)

            print("{0} (x/t <= {1:.3e})".format(self.left, self.sl))
            print("{0} ({1:.3e} < x/t <= {2:.3e})".format(self.wlstar, self.sl, self.vstar))
        else:
            self.densl = self.left.dens * math.pow(self.pstar / self.left.pres, 1.0 / self.gamma)
            al = self.left.sonic * math.pow(self.pstar / self.left.pres, (self.gamma - 1.0) / (2.0 * self.gamma))
            self.shl = self.left.vel - self.left.sonic
            self.stl = self.vstar - al
            self.wlstar = Wave(self.vstar, self.densl, self.pstar, self.gamma)
            print("{0} (x/t <= {1:.3e})".format(self.left, self.shl))
            print("rarefaction wave                                                ({0:.3e} < x/t <= {1:.3e})".format(self.shl, self.stl))
            print("{0} ({1:.3e} < x/t <= {2:.3e})".format(self.wlstar, self.stl, self.vstar))


        if self.right_shock:
            self.densr = self.right.dens * ((self.gamma - 1.0) * self.right.pres + (self.gamma + 1.0)
                    * self.pstar) / ((self.gamma + 1.0) * self.right.pres + (self.gamma - 1.0) * self.pstar)
            self.sr = self.right.vel + self.right.sonic * math.sqrt((self.gamma + 1.0) / (2.0 * self.gamma)
                    * self.pstar / self.right.pres + (self.gamma - 1.0) / (2.0 * self.gamma))
            self.wrstar = Wave(self.vstar, self.densr, self.pstar, self.gamma)
            print("{0} ({1:.3e} < x/t <= {2:.3e})".format(self.wrstar, self.vstar, self.sr))
            print("{0} ({1:.3e} < x/t)".format(self.right, self.sr))
        else:
            self.densr = self.right.dens * math.pow(self.pstar / self.right.pres, 1.0 / self.gamma)
            ar = self.right.sonic * math.pow(self.pstar / self.right.pres, (self.gamma - 1.0) / (2.0 * self.gamma))
            self.shr = self.right.vel + self.right.sonic
            self.str = self.vstar + ar
            self.wrstar = Wave(self.vstar, self.densr, self.pstar, self.gamma)
            print("{0} ({1:.3e} < x/t <= {2:.3e})".format(self.wrstar, self.vstar, self.str))
            print("rarefaction wave                                                ({0:.3e} < x/t <= {1:.3e})".format(self.str, self.shr))
            print("{0} ({1:.3e} < x/t)".format(self.right, self.shr))

    def output(self, num, time, filename):
        dx = 1.0 / num
        w = []
        for i in range(num):
            x = -0.5 + (0.5 + i) * dx
            x_t = x / time

            if self.left_shock:
                if x_t <= self.sl:
                    w.append(self.left)
                    continue
                elif self.sl < x_t <= self.vstar:
                    w.append(self.wlstar)
                    continue
            else:
                if x_t <= self.shl:
                    w.append(self.left)
                    continue
                elif self.shl < x_t <= self.stl:
                    tmp = (2.0 / (self.gamma + 1.0) + (self.gamma - 1.0) / ((self.gamma + 1.0) * self.left.sonic) * (self.left.vel - x_t))
                    wfan = Wave(
                            2.0 / (self.gamma + 1.0) * (self.left.sonic + (self.gamma - 1.0) * 0.5 * self.left.vel + x_t),
                            self.left.dens * math.pow(tmp, 2.0 / (self.gamma - 1.0)),
                            self.left.pres * math.pow(tmp, 2.0 * self.gamma / (self.gamma - 1.0)),
                            self.gamma)
                    w.append(wfan)
                    continue
                elif self.stl < x_t <= self.vstar:
                    w.append(self.wlstar)
                    continue

            if self.right_shock:
                if x_t >= self.sr:
                    w.append(self.right)
                    continue
                elif self.sr > x_t >= self.vstar:
                    w.append(self.wrstar)
                    continue
            else:
                if x_t >= self.shr:
                    w.append(self.right)
                    continue
                elif self.shr > x_t >= self.str:
                    tmp = (2.0 / (self.gamma + 1.0) - (self.gamma - 1.0) / ((self.gamma + 1.0) * self.right.sonic) * (self.right.vel - x_t))
                    wfan = Wave(
                            2.0 / (self.gamma + 1.0) * (-self.right.sonic + (self.gamma - 1.0) * 0.5 * self.right.vel + x_t),
                            self.right.dens * math.pow(tmp, 2.0 / (self.gamma - 1.0)),
                            self.right.pres * math.pow(tmp, 2.0 * self.gamma / (self.gamma - 1.0)),
                            self.gamma)
                    w.append(wfan)
                    continue
                elif self.str > x_t >= self.vstar:
                    w.append(self.wrstar)
                    continue

        with open(filename, "w") as f:
            f.write("# time {0}\n".format(time))
            for i, v in enumerate(w):
                x = (i + 0.5) * dx - 0.5
                f.write("{0} {1} {2} {3} {4}\n".format(x, v.vel, v.dens, v.pres, v.pres / ((self.gamma - 1.0) * v.dens)))

        print("output {0}".format(filename))

def main():
    args = sys.argv
    if len(args) == 2:
        with open(args[1], "r") as f:
            json_data = json.load(f)
            l = json_data["left"]
            r = json_data["right"]
            gamma = json_data["gamma"]
            num = json_data["number"]
            time = json_data["time"]
            fname = json_data["fileName"]

            left  = Wave(l["velocity"], l["density"], l["pressure"], gamma)
            right = Wave(r["velocity"], r["density"], r["pressure"], gamma)

        rs = RiemannSolver(left, right, gamma)
        rs.solve()
        rs.output(num, time, fname)

    else:
        print("how to use")
        print("riemann_solver.py <json file>")

if __name__ == "__main__":
    main()
