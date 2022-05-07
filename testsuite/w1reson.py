#! /usr/bin/env python

#############################################################################
#
# W1 waveguide.
#
#############################################################################

from camfr import *

import unittest, eps

class w1reson(unittest.TestCase):
    def create_w1_period(self,period,rad,periods,topbuf,matcore,matclad):
        tussen=sqrt(3)/2 *period
        tussen2 = sqrt(3)*period
        numodd = int(periods/2)
        numeven = int((periods+1)/2)
        propstart,propend,propprec = 0,2*rad+period/2,period/5
        transstart,transend,transprec = 0,tussen*periods+rad+topbuf,period/5
        gi = Geometry(matcore)
        for i in range(0,numodd,1):
                gi+= Circle(Point(rad,tussen+tussen2*i),rad,matclad)
                gi+= Circle(Point(rad+period,tussen+tussen2*i),rad,matclad)
        for i in range(1,numeven+1,1):
                gi += Circle(Point(rad+period/2.0,tussen2*i),rad,matclad)
        return gi.to_expression(rad,rad+period,propprec,transstart,transend,transprec)

    def create_w1_1def_period(self,period,rad,defrad,periods,topbuf,matcore,matclad):
        tussen=sqrt(3)/2 *period
        tussen2 = sqrt(3)*period
        numodd = int(periods/2)
        numeven = int((periods+1)/2)
        propstart,propend,propprec = 0,2*rad+period/2,period/30
        transstart,transend,transprec = 0,tussen*periods+rad+topbuf,period/30
        gi = Geometry(matcore)
        for i in range(0,numodd,1):
                gi+= Circle(Point(rad,tussen+tussen2*i),rad,matclad)
                gi+= Circle(Point(rad+period,tussen+tussen2*i),rad,matclad)
        for i in range(1,numeven+1,1):
                gi += Circle(Point(rad+period/2.0,tussen2*i),rad,matclad)
        gi += Circle(Point(rad+period/2.0,0.0),defrad,matclad)
        return gi.to_expression(rad,rad+period,propprec,transstart,transend,transprec)

    def testw1reson(self):
        
        """W1 resonator"""

        print("")
        print("Running W1 resonator...")
        
        set_lambda(1.526)

        set_N(15)
        set_mode_surplus(5)
        set_solver(series)
        set_orthogonal(0)
        set_polarisation(TM)

        SOI = Material(2.7)
        Air = Material(1.0)
        a = 0.47
        r = 0.3*a
        periods = 2
        topbuf = -r+a*0.5*sqrt(3)-0.25

        set_lower_wall(slab_E_wall)
        set_upper_PML(0)
        set_lower_PML(0)

        w1 = self.create_w1_period(a,r,periods,topbuf,SOI,Air)
        w1d = self.create_w1_1def_period(a,r,0.95*r,periods,topbuf,SOI,Air)

        wgb = BlochStack(w1)
        wgb.calc()
        wg = Stack(wgb(0)+w1*10+w1d+w1*30+w1d+w1*10+wgb(0))
        wg.calc()

        t = wg.T12(0,0)
        r = wg.R12(0,0)
    
        ta = abs(t)**2
        ra = abs(r)**2
        tp = atan(t.imag/t.real)
        rp = atan(r.imag/r.real)

        ta_OK =  0.235682188665
        tp_OK = -0.46442242995
        ra_OK =  0.764317811335
        rp_OK =  0.291072912404

        print(ta, "expected", ta_OK)
        print(tp, "expected", tp_OK)
        print(ra, "expected", ra_OK)
        print(rp, "expected", rp_OK)

        ta_pass = abs((ta-ta_OK)/ta_OK) < eps.testing_eps
        tp_pass = abs((tp-tp_OK)/tp_OK) < eps.testing_eps
        ra_pass = abs((ra-ra_OK)/ra_OK) < eps.testing_eps
        rp_pass = abs((rp-rp_OK)/rp_OK) < eps.testing_eps
    
        free_tmps()

        self.assertTrue(ta_pass and tp_pass and ra_pass and rp_pass)

        set_solver(track)
        set_orthogonal(1)
        set_polarisation(TE)

suite = unittest.makeSuite(w1reson, 'test')        

if __name__ == "__main__":
    unittest.main()
    




