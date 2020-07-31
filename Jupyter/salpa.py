#!/usr/bin/python3

import numpy as np

ASYM_NOT_CHI2 = True
PREMATURE = True
THIRDORDER = True
TOOPOORCNT = 5

class State:
    OK = 0
    PEGGING = 1
    PEGGED = 2
    TOOPOOR = 3
    DEPEGGING = 4
    FORCEPEG = 5
    BLANKDEPEG = 6

class Salpa:
    def update_X012(self):
        y_new = self.source[self.t_stream + self.tau]
        y_old = self.source[self.t_stream - self.tau - 1]
        self.X0 += y_new - y_old
        self.X1 += self.tau_plus_1*y_new - self.minus_tau*y_old - self.X0
        self.X2 += self.tau_plus_1_squared*y_new \
                   - self.minus_tau_squared*y_old - self.X0 - 2*self.X1

    def calc_alpha0(self):
        self.alpha0 = (self.T4*self.X0 - self.T2*self.X2) \
                      / (self.T0*self.T4-self.T2*self.T2)

    def __init__(self,
                 thresh=np.inf, tau=30,
                 t_blankdepeg=5, t_ahead=5, t_chi2=15,
                 rails=[-np.inf, np.inf]):
        self.y_threshold = thresh
        self.tau = tau
        self.t_blankdepeg = t_blankdepeg
        self.t_ahead = t_ahead
        self.t_chi2 = t_chi2
        self.rail1 = rails[0]
        self.rail2 = rails[1]

    def reset(self, t_start):
        self.t_peg = t_start
        self.t_stream = t_start
        self.t0 = None
        self.state = State.PEGGED

    def ispegged(self, v):
        return v<=self.rail1 or v>=self.rail2

    def init_T(self):
        if ASYM_NOT_CHI2:
            self.my_thresh = 3.92 * self.t_chi2 * self.y_threshold**2
        else:
            self.my_thresh = (self.t_chi2 - 4) * self.y_threshold**2

        self.tau_plus_1 = self.tau + 1
        self.tau_plus_1_squared = self.tau_plus_1**2
        self.tau_plus_1_cubed = self.tau_plus_1**3
        self.minus_tau = -self.tau
        self.minus_tau_squared = self.minus_tau**2
        self.minus_tau_cubed = self.minus_tau**3
        self.T0 = 0
        self.T2 = 0
        self.T4 = 0
        self.T6 = 0
        for t in range(-self.tau, self.tau+1):
            self.T0 += 1
            self.T2 += t**2
            self.T4 += t**4
            self.T6 += t**6

    def apply(self, source):
        self.source = source
        self.dest = np.zeros(source.shape, source.dtype)
        self.t_end = len(self.source)
        self.init_T()
        self.reset(0)
        self.process(self.t_end)
        res = self.dest
        self.source = None
        self.dest = None
        return res

    def process(self, t_limit):
        self.state = self.statemachine(t_limit, self.state)
        return self.t_stream

    def forcepeg(self, t_from, t_to):
        self.state = self.statemachine(t_from - self.tau, self.state)
        if self.state==State.OK:
            self.t0 = self.t_stream - 1
            self.calc_X3()
            self.calc_alpha0123()
            self.statemachine(t_from, State.PEGGING)
        self.t0 = t_to
        self.state = self.statemachine(t_to, State.FORCEPEG)
        return self.t_stream

    def state_pegged(self):
        if self.ispegged(self.source[self.t_stream]):
            self.dest[self.t_stream] = np.nan
            self.t_stream += 1
            return State.PEGGED
        for dt in range(1, 2*self.tau+1):
            if self.t_stream + dt >= self.t_end \
               or self.ispegged(self.source[self.t_stream + dt]):
                self.t0 = self.t_stream + dt
                return State.FORCEPEG
        self.t0 = self.t_stream + self.tau
        self.calc_X012()
        self.calc_X3()
        self.calc_alpha0123()
        self.toopoorcnt = TOOPOORCNT
        return State.TOOPOOR

    def state_toopoor(self):
        if ASYM_NOT_CHI2:
            asym = 0
            sig = 0
            for i in range(self.t_chi2):
                t_i = self.t_stream + i
                if t_i > self.t_end:
                    self.t0 = t_i
                    return State.FORCEPEG
                dt = t_i - self.t0
                dy = self.alpha0 + self.alpha1*dt + self.alpha2*dt**2 \
                     + self.alpha3*dt**3 - self.source[t_i]
                asym += dy
                sig += dy*dy
            asym *= asym
            if asym < self.my_thresh:
                self.toopoorcnt -= 1
                if self.toopoorcnt<=0 and asym<self.my_thresh/3.92:
                    if PREMATURE:
                        dt = self.t_stream - self.t0
                        self.negv = self.source[self.t_stream] \
                                    < self.alpha0 + self.alpha1*dt \
                                    + self.alpha2*dt**2 + self.alpha3*dt**3
                    self.calc_X012()
                    self.calc_X3()
                    return State.BLANKDEPEG
            else:
                self.toopoorcnt = TOOPOORCNT
        else:
            chi2 = 0
            for i in range(self.t_chi2):
                t_i = self.t_stream + self.t_blankdepeg + i
                if t_i > self.t_end:
                    self.t0 = t_i
                    return State.FORCEPEG
                dt = t_i - self.t0
                dy = self.alpha0 + self.alpha1*dt \
                     + self.alpha2*dt**2 + self.alpha3*dt**3\
                     - self.source[t_i]
                chi2 += dy*dy
            if chi2 < self.my_thresh:
                if PREMATURE:
                    dt = self.t_stream - self.t0
                    self.negv = self.source[self.t_stream] \
                                < self.alpha0 + self.alpha1*dt \
                                + self.alpha2*dt**2 + self.alpha3*dt**3
                return State.BLANKDEPEG
        self.dest[self.t_stream] = np.nan
        self.t_stream += 1
        self.t0 += 1
        t1 = self.t0 + self.tau
        if t1 >= self.t_end or self.ispegged(self.source[t1]):
            self.t0 += self.tau
            return State.FORCEPEG
        self.update_X0123() # Pointless?
        self.calc_X012()
        self.calc_X3()
        self.calc_alpha0123()
        return State.TOOPOOR

    def state_forcepeg(self):
        if self.t_stream >= self.t0:
            return State.PEGGED
        self.dest[self.t_stream] = np.nan
        self.t_stream += 1
        return State.FORCEPEG

    def state_blankdepeg(self):
        if self.t_stream >= self.t0 - self.tau + self.t_blankdepeg:
            return State.DEPEGGING
        if PREMATURE:
            dt = self.t_stream - self.t0
            y = self.source[self.t_stream] - (self.alpha0
                                              + self.alpha1*dt
                                              + self.alpha2*dt**2
                                              + self.alpha3*dt**3)
            #print('blankdepeg', self.t0, self.t_stream, y)
            if (y<0) != self.negv:
                self.dest[self.t_stream] = y;
                self.t_stream += 1
                return State.DEPEGGING
        self.dest[self.t_stream] = np.nan
        self.t_stream += 1
        return State.BLANKDEPEG

    def state_depegging(self):
        if self.t_stream==self.t0:
            return State.OK
        dt = self.t_stream - self.t0
        self.dest[self.t_stream] = self.source[self.t_stream] \
                                   - (self.alpha0 + self.alpha1*dt
                                      + self.alpha2*dt**2 + self.alpha3*dt**3)
        self.t_stream += 1
        return State.DEPEGGING

    def state_pegging(self):
        if self.t_stream >= self.t0 + self.tau:
            self.t_peg = self.t_stream
            return State.PEGGED
        dt = self.t_stream - self.t0
        self.dest[self.t_stream] = self.source[self.t_stream] \
                                   - (self.alpha0 + self.alpha1*dt
                                      + self.alpha2*dt**2 + self.alpha3*dt**3)
        self.t_stream += 1
        return State.PEGGING
        
    def state_ok(self):
        self.calc_alpha0()
        self.dest[self.t_stream] = self.source[self.t_stream] - self.alpha0
        self.t_stream += 1
        t1 = self.t_stream + self.tau + self.t_ahead
        if t1 >= self.t_end or self.ispegged(self.source[t1]):
            self.t0 = self.t_stream - 1
            self.calc_X3()
            self.calc_alpha0123()
            return State.PEGGING
        self.update_X012()
        return State.OK
    
    statemap = { State.OK: state_ok,
                 State.PEGGED: state_pegged,
                 State.PEGGING: state_pegging,
                 State.TOOPOOR: state_toopoor,
                 State.DEPEGGING: state_depegging,
                 State.FORCEPEG: state_forcepeg,
                 State.BLANKDEPEG: state_blankdepeg
                 }
    
    def statemachine(self, t_limit, s):
        while self.t_stream < t_limit:
            #print(self.t_stream, self.t0, s)
            s = self.statemap[s](self)
        return s

    def calc_X012(self):
        self.X0 = 0
        self.X1 = 0
        self.X2 = 0
        for t in range(-self.tau, self.tau+1):
            y = self.source[self.t0 + t]
            self.X0 += y
            self.X1 += t * y
            self.X2 += t**2 * y

    def calc_X3(self):
        self.X3 = 0
        for t in range(-self.tau, self.tau+1):
            y = self.source[self.t0 + t]
            self.X3 += t**3 * y

    def update_X0123(self):
        y_new = self.source[self.t0 + self.tau]
        y_old = self.source[self.t0 - self.tau - 1]
        self.X0 += y_new - y_old
        self.X1 += self.tau_plus_1*y_new - self.minus_tau*y_old - self.X0
        self.X2 += self.tau_plus_1_squared*y_new - self.minus_tau_squared*y_old \
                   - self.X0 - 2*self.X1
        self.X3 += self.tau_plus_1_cubed*y_new - self.minus_tau_cubed*y_old \
                   - self.X0 - 3*self.X1 - 3*self.X2

    def calc_alpha0123(self):
        fact02 = 1./(self.T0*self.T4 - self.T2*self.T2)
        self.alpha0 = fact02*(self.T4*self.X0 - self.T2*self.X2)
        self.alpha2 = fact02*(self.T0*self.X2 - self.T2*self.X0)
        if THIRDORDER:
            fact13 = 1./(self.T2*self.T6-self.T4*self.T4)
            self.alpha1 = fact13*(self.T6*self.X1 - self.T4*self.X3)
            self.alpha3 = fact13*(self.T2*self.X3 - self.T4*self.X1)
        else:
            self.alpha1 = self.X1/self.T2
            self.alpha3 = 0
        #self.report()
        
    def report(self):
        print(f'state={self.state} t_stream={self.t_stream} t0={self.t0} y[t]={self.source[self.t_stream]:.4} alpha={self.alpha0:.4} {self.alpha1:.4} {self.alpha2:.4} {self.alpha3:.4} X={self.X0:.4} {self.X1:.4} {self.X2:.4} {self.X3:.4}')
        
