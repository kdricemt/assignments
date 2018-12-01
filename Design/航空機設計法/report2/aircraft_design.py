%matplotlib inline

import math as m
import numpy as np
import matplotlib.pyplot as plt

#this sizing method is only appliable for jet-engine aircrafts
class Sizing:
    rho_0 = 0.125 #kgs^2/m^4
    a_0 = 340.3 #[m/s]
    #unit conversion
    m2knot_0 = 1.943
    knot2feets = 1.687
    kgm4toibft4 = 2.204/(3.280**4)
    m2ft = 3.280

    def __init__(self):
    # parameters required
        # passengers
        self.passengers = 420
        self.pilots = 2
        # flight parameter
        self.fr = 7500 #flight range
        self.fr_alt = 200 #for alternative airport
        self.e_ltr = 0.75 # [hr]
        self.c_j = 0.45 #specific fuel consumption
        self.a = 0.4736 # log10(W_TO) = A + Blog10(W_E)
        self.b = 0.9656
        self.stofl = 10000 #take-off field length
        self.sfl = 7000 # Far landing field length
        self.wlbywto_sfl = 0.8 # WL/W_TO at sealevel
        self.lapserate_cruise = 0.167 #Lp cr(p100)
        self.altitude = 38000 #[ft]
        self.rho_cruise = 0.272 * Sizing.rho_0 #density at cruising altitude
        self.a_cruise = 0.867 * Sizing.a_0 #a at cruising altitude
        self.m_cruise = 0.8 #mach number cruise
        self.v = self.m_cruise * self.a_cruise * Sizing.m2knot_0 #cruise speed [m/s]

    #parameters to set
        self.economy_ratio = 0.84 #the ratio of passengers in economy class
        self.buisness_ratio = 0.10
        self.lbyd_guess = 20 # L/D used in W_TO calculation
        self.lbyd_ltr = 23

        self.ar = 9.0 #aspect ratio
        self.swetbysref = 3 #Swet/S
        self.war = self.ar / self.swetbysref #wetted aspect ratio
        self.lbydmax = 23 # p81 graph
        self.n_engines = 3 #number of engines

        #plot parameter
        self.wbysto_min = 30 #input array range
        self.wbysto_max = 220
        self.xmax = 250
        self.ymax = 0.6
        #array for inputs
        self.wbysto_array_input = np.linspace(self.wbysto_min,self.wbysto_max,self.wbysto_max - self.wbysto_min)
        self.clmax_l_array_input= [1.9,2.2,2.6,3.0]
        self.clmax_to_array_input = [1.7,2.0,2.4]
        self.clmax_to_array_input_lift = 2.4

    #parameters calculated
        self.economy = 0 #passengers in economy class
        self.buisness = 0
        self.first = 0
        self.ca = 0 #cabin attendant
        self.crew = 0

        self.w5byw4 = 0
        self.w7byw6 = 0
        self.w8byw7 = 0
        self.m_ff = 0 #mission fuel function
        self.w_to = 0 #maximum take-off weight
        self.w_crew = 0 #crew
        self.w_pl = 0 #payloads
        self.w_etent = 0
        self.w_oetent = 0
        self.w_e = 0
        self.w_error = 0
        self.lbyd_cal = 0

        self.cd0 = 0
        self.cd0_to = 0 #Take-off
        self.cd0_toleg = 0 #Take-off legs out
        self.cd0_la = 0 #Landing
        self.cd0_laleg = 0

        self.epiar = 0
        self.e_to = 0
        self.e_la = 0

        self.vsl = 0

        #design point
        self.wbysto_designpoint = 0
        self.tbywto_designpoint = 0
        self.thrust = 0 # Total thrust
        self.wing_area = 0

    #display functions
    def display_requirements(self):
        print('======================================================')
        print('Required Parameters')
        print('  Passengers               ',self.passengers)
        print('  Flight Range             ',self.fr,'[nm]')
        print('    (for alternate airport)',self.fr_alt,'[nm]','(',self.e_ltr*60,'min)')
        print('  Cruise Speed             ','{:.2f}'.format(self.v),'[m/s] ','M',self.m_cruise)
        print('  Altitude                 ',self.altitude,'[ft]')
        print('  C_j                      ',self.c_j)
        print('  Take-off field length    ',self.stofl,'[ft]')
        print('  FAR landing field length ',self.sfl,'[ft]  at (W_L/W_TO)=',self.wlbywto_sfl)
        print('  Lapse rate_cruise        ',self.lapserate_cruise)
        print('  log10(W_TO) =',self.a,'+',self.b,'log10(W_E)')

    def display_setparameter(self):
        print('======================================================')
        print('Setting Parameters')
        print('  (L/D) Guess              ',self.lbyd_guess)
        print('  (L/D) ltr                ',self.lbyd_ltr)
        print('')

        print('  AR                       ',self.ar)
        print('  Wetted Aspect Ratio      ',self.war)
        print('  (L/D) Max                ',self.lbydmax)
        print('  Number of Engines        ',self.n_engines)

    def display_results(self):
        print('======================================================')
        print('Results')
        print('  <W_crew and W_PL>---------------------------------p75')
        print('     Passengers            ',self.passengers)
        print('         Economy           ',self.economy)
        print('         Buisness          ',self.buisness)
        print('         First             ',self.first)
        print('     Crew                  ',self.crew)
        print('          Cabin Attendant  ',self.ca)
        print('          Pilots           ',self.pilots)
        print('     W_crew                ','{:.0f}'.format(self.w_crew),'[Ibs]')
        print('     W_pl                  ','{:.0f}'.format(self.w_pl),'[Ibs]')
        print('')

        print('  <Maximum Take-off Weight>------------------------p65')
        print('     W5/W4                 ','{:3f}'.format(self.w5byw4))
        print('     W7/W6                 ','{:3f}'.format(self.w7byw6))
        print('     W8/W7                 ','{:3f}'.format(self.w8byw7))
        print('     M_ff                  ','{:.2f}'.format(self.m_ff),'[Ibs]')
        print('     W_TO                  ','{:.0f}'.format(self.w_to),'[Ibs]')
        print('     W_OEtent              ','{:.0f}'.format(self.w_oetent),'[Ibs]')
        print('     W_Etent               ','{:.0f}'.format(self.w_etent),'[Ibs]')
        print('     W_E                   ', '{:.0f}'.format(self.w_e),'[Ibs]')
        print('     W_error               ','{:.6f}'.format(self.w_error))
        print('')

        print('  <L/D and Drug>-----------------------------------p81')
        print('     (L/D)max              ','{:.1f}'.format(self.lbydmax))
        print('     (L/D)cruise           ','{:.1f}'.format(self.lbyd_cal))
        print('     Clean                 ','{:.4f}'.format(self.cd0),'+' \
        ,'{:.4f}'.format(self.epiar),'CL^2' )
        print('     Take-off Gear-up      ','{:.4f}'.format(self.cd0_to),'+' \
        ,'{:.4f}'.format(self.epiar_to),'CL^2' )
        print('     Take-off Gear-down    ','{:.4f}'.format(self.cd0_toleg),'+', \
        '{:.4f}'.format(self.epiar_to),'CL^2' )
        print('     Landing Gear-up       ','{:.4f}'.format(self.cd0_la),'+', \
        '{:.4f}'.format(self.epiar_la),'CL^2' )
        print('     Landing Gear-down     ','{:.4f}'.format(self.cd0_laleg),'+', \
        '{:.4f}'.format(self.epiar_la),'CL^2' )

    #calculation functions
    def calculation(self):
        self._cal_payloadcrew()
        self._cal_mff()
        self._cal_w_to()
        self._cal_ld_drug()

    #calculate maximum lift-off weight
    def _cal_payloadcrew(self):
      self.economy = int(self.passengers*self.economy_ratio)
      self.buisness = int(self.passengers*self.buisness_ratio)
      self.first = self.passengers - (self.economy + self.buisness)
      self.ca = int(self.economy/30) + int(self.buisness/20) + int(self.first/10)
      self.crew = self.ca + self.pilots

      self.w_crew = self.crew*(175+30)
      self.w_pl = self.passengers*175 + self.economy*44 + (self.buisness+self.first)*66

    def _cal_mff(self):
        # product of weight ratio aready known
        wr_known = 0.990*0.990*0.995*0.980*0.990*0.992
        # for alternative airport
        v_alt = 300
        lbyd_alt = self.lbyd_guess
        cj_alt = self.c_j
        #during loiter
        cj_ltr = 0.4

        self.w5byw4 = m.exp(-self.fr / (self.v * self.lbyd_guess / self.c_j))
        self.w7byw6 = m.exp(-self.fr_alt / (self.v * lbyd_alt / cj_alt))
        self.w8byw7 = m.exp(-self.e_ltr / (self.lbyd_ltr / cj_ltr))
        mff = wr_known * self.w5byw4 * self.w7byw6 * self.w8byw7
        self.m_ff = mff

    #find W_TO that makes W_E - W_Etent smallest
    def _cal_w_to(self):
        counter = 0
        start = 500
        end = 1500
        w_oetent_guess = np.zeros(end-start)
        w_etent_guess = np.zeros(end-start)
        w_e_guess = np.zeros(end-start)
        w_diff = np.zeros(end-start)
        wto_guess = np.zeros(end-start)
        w_fused_ratio = 1-self.m_ff #W_fused/W_to

        for wto1000 in range(start,end):
            wto_guess[counter] = wto1000 * 1000
            w_oetent_guess[counter] = (1-w_fused_ratio)*wto_guess[counter] - self.w_pl
            w_etent_guess[counter] = w_oetent_guess[counter] - self.w_crew
            w_e_guess[counter] = m.pow(10,(m.log10(wto_guess[counter])-self.a)/self.b)
            w_diff[counter] = abs(w_etent_guess[counter] - w_e_guess[counter])
            counter = counter + 1

        min_index = np.argmin(w_diff)
        self.w_to = wto_guess[min_index]
        self.w_oetent = w_oetent_guess[min_index]
        self.w_etent = w_etent_guess[min_index]
        self.w_e = w_e_guess[min_index]
        self.w_error = abs(w_diff[min_index])/self.w_e

    def _epiar(self,e):
        epiar = 1/(e * m.pi * self.ar)
        return epiar

    def _cal_ld_drug(self):
        cfe = 0.003 #equivlalent skin friction coefficient
        e = 0.80

        self.lbyd_cal = self.lbydmax * 0.886

        self.cd0 = cfe * self.swetbysref
        self.cd0_to = self.cd0 + 0.015
        self.cd0_toleg = self.cd0 + 0.015 + 0.020
        self.cd0_la = self.cd0 + 0.065
        self.cd0_laleg = self.cd0 + 0.065 + 0.020

        self.epiar = self._epiar(e)
        self.epiar_to = self._epiar(e + 0.75)
        self.epiar_la = self._epiar(e + 0.70)

    #sizing for take-off
    def _tbyw_to_takeoff(self,clmax_to,wbys_to):
        tbyw = 40.3*wbys_to/(self.stofl*clmax_to)
        if wbys_to == self.wbysto_max:
            print('    CL_max_TO',clmax_to,' (T/W)TO = ', '{:.5f}'.format(40.3/(self.stofl*clmax_to)), '(W/S)TO')
        return tbyw

    #sizing for landing
    def _cal_vsl(self):
        self.vsl = m.sqrt(self.sfl/0.29)/1.3 #[knot]

    def _wbys_to_landing(self,cl_max_l):
        self._cal_vsl()
        #convert vsl to ft/s
        wbys_to = 0.5*(Sizing.rho_0 * Sizing.kgm4toibft4)*((self.vsl * Sizing.knot2feets)**2) * cl_max_l / self.wlbywto_sfl
        print('    CL_max_l ',cl_max_l,' (W/S)TO = ','{:.1f}'.format(wbys_to),'[Ib/ft^2]')
        if cl_max_l == self.clmax_l_array_input[3]:
            print('    V_SL               ','{:.0f}'.format(self.vsl),'[kt]')
        return wbys_to

    #sizing for lift
    def _tbyw_to_min_lift(self,cl_max_to):
        gamma = 0.024 + (self.n_engines-2)*0.003
        vbyvsto = 1.2 # V = 1.2 V_STO
        cl_lift = cl_max_to/(vbyvsto**2)
        cd_lift = self.cd0_to + self.epiar_to*(cl_lift**2) #liftoff no-leg
        lbyd_lift = cl_lift/cd_lift
        print('    CD                 ','{:.3f}'.format(cd_lift))
        print('    L/D                ','{:.1f}'.format(lbyd_lift))
        tbyw_to_min = (self.n_engines/(self.n_engines-1))*((1/lbyd_lift)+gamma)
        print('    (T/W)_TO           ','{:.3f}'.format(tbyw_to_min))
        tbyw_to_min = tbyw_to_min/0.8
        print('    (T/W)_TO (/0.8)    ','{:.3f}'.format(tbyw_to_min))
        return tbyw_to_min

    #sizing for cruise
    def _tbyw_to_cruise(self,wbys_to):
        delta_cd0 = 0.003
        q = 0.5 * (self.rho_cruise * Sizing.kgm4toibft4) * ((self.v * Sizing.knot2feets)**2)
        w_crbyw_to = 0.956
        wbys_cr = 0.956 * wbys_to
        tbyw_cr = ((self.cd0 + delta_cd0) * q / wbys_cr) + (wbys_cr * self.epiar/ q)
        tbywto = tbyw_cr * w_crbyw_to / self.lapserate_cruise
        if wbys_to == self.wbysto_max:
            print('    q                  ','{:.0f}'.format(q), '[Ib/ft^2]')
            print('    (T/W)_TO=', '{:.2f}'.format(w_crbyw_to / self.lapserate_cruise),'(', \
            '{:.2f}'.format((self.cd0 + delta_cd0) * q / w_crbyw_to)  ,'/(W/S)_TO + (W/S)_TO/', \
            '{:.0f}'.format((1/w_crbyw_to) * q * (1/self.epiar)) ,')')
        return tbywto

    #plot graphs
    def plot(self):
        #vectorize
        v_tbyw_to_takeoff = np.vectorize(self._tbyw_to_takeoff) #sizing for take-off
        v_wbys_to_landing = np.vectorize(self._wbys_to_landing) #sizing for landing
        v_tbyw_to_cruise = np.vectorize(self._tbyw_to_cruise)

        #array for inputs
        #wbysto_array_input = np.linspace(self.wbysto_min,self.wbysto_max,self.wbysto_max - self.wbysto_min)
        #clmax_l_array_input= [1.8,2.2,2.6,3.0]
        #clmax_to_array_input = [1.6,2.0,2.4]
        #clmax_to_array_input_lift = 2.4

        #array for results
        #take off
        print('')
        print('  <Sizing for Take-off>---------------------------p85')
        tbyw_to_takeoff_array0 = v_tbyw_to_takeoff(self.clmax_to_array_input[0],self.wbysto_array_input)
        tbyw_to_takeoff_array1 = v_tbyw_to_takeoff(self.clmax_to_array_input[1],self.wbysto_array_input)
        tbyw_to_takeoff_array2 = v_tbyw_to_takeoff(self.clmax_to_array_input[2],self.wbysto_array_input)
        #landing
        print('')
        print('  <Sizing for Landing>----------------------------p89')
        wbys_to_landing_array = v_wbys_to_landing(self.clmax_l_array_input) #parrael to y
        #lift
        print('')
        print('  <Sizing for Lift>-------------------------------p91')
        tbyw_to_min_lift = self._tbyw_to_min_lift(self.clmax_to_array_input_lift) #parrarel to x
        #cruise
        print('')
        print('  <Sizing for Cruise>-----------------------------p98')
        tbyw_to_cruise_array = v_tbyw_to_cruise(self.wbysto_array_input)
        print('=======================================================')
        print('Sizing Result (at design point P)')
        print('')

        #print result arrays
        #print('tbyw_to_takeoff_array0:',tbyw_to_takeoff_array0)
        #print('wbys_to_landing_array:',wbys_to_landing_array)
        #print('tbyw_to_min_lift_array:',tbyw_to_min_lift_array)
        #print('tbyw_to_cruise_array:',tbyw_to_cruise_array)

        #plot
        fig = plt.figure()
        plt.xlim(0,self.xmax)
        plt.ylim(0,self.ymax)
        ax1 = fig.add_subplot(111)
        ax2 = fig.add_subplot(111)
        ax3 = fig.add_subplot(111)
        ax4 = fig.add_subplot(111)
        ax5 = fig.add_subplot(111)

        ax1.plot(self.wbysto_array_input,tbyw_to_takeoff_array0,color='black',linewidth=1.0)
        ax2.plot(self.wbysto_array_input,tbyw_to_takeoff_array1,color='black',linewidth=1.0)
        ax3.plot(self.wbysto_array_input,tbyw_to_takeoff_array2,color='black',linewidth=1.0)
        ax4.plot(self.wbysto_array_input,tbyw_to_cruise_array,color='black',linewidth=1.0)
        ax5.scatter(self.wbysto_designpoint,self.tbywto_designpoint,color='black')

        plt.vlines([wbys_to_landing_array],0,self.ymax-0.1,color="black",linewidth=0.5)
        plt.hlines([tbyw_to_min_lift],0,self.xmax,color='black',linewidth=0.5)

        ax1.set_xlabel("$ (W/S)_{TO} \ [Ib/ft^2]$")
        ax1.set_ylabel("$ (T/W)_{TO} $")

        plt.text(wbys_to_landing_array[0]-40,0.52, r'$C_{LmaxL} =$', fontsize = 9)
        plt.text(wbys_to_landing_array[0]-5,0.52,self.clmax_l_array_input[0], fontsize = 9)
        plt.text(wbys_to_landing_array[1]-5,0.52,self.clmax_l_array_input[1], fontsize = 9)
        plt.text(wbys_to_landing_array[2]-5,0.52,self.clmax_l_array_input[2], fontsize = 9)
        plt.text(wbys_to_landing_array[3]-5,0.52,self.clmax_l_array_input[3], fontsize = 9)
        plt.text(5,0.30, 'Cruise Speed', fontsize=9)
        plt.text(self.wbysto_max,tbyw_to_takeoff_array0[-1],r'$C_{LmaxTO}$', fontsize = 9)
        str0 = str(self.clmax_to_array_input[0])
        str1 = str(self.clmax_to_array_input[1])
        str2 = str(self.clmax_to_array_input[2])
        plt.text(self.wbysto_max,tbyw_to_takeoff_array0[-1]-0.03,'='+str0, fontsize = 9)
        plt.text(self.wbysto_max,tbyw_to_takeoff_array1[-1],'='+str1, fontsize = 9)
        plt.text(self.wbysto_max,tbyw_to_takeoff_array2[-1],'='+str2, fontsize = 9)
        plt.text(5,tbyw_to_min_lift+0.02,'FAR25.121(OEI)', fontsize=9)
        plt.text(self.wbysto_designpoint+5,self.tbywto_designpoint,'P', fontsize=11)
        filename = "sizing.png"
        fig.savefig(filename)

        print('    (W/S)_TO             ','{:.0f}'.format(self.wbysto_designpoint))
        print('    (T/W)_TO             ','{:.4f}'.format(self.tbywto_designpoint))

        print('    S                    ','{:.1f}'.format(self.wing_area),' [ft^2]')
        print('    T                    ','{:.0f}'.format(self.thrust),' [Ib]')

    def  cal_s_t(self,clmax_l_designpoint,clmax_to_designpoint):
        self.wbysto_designpoint = self._wbys_to_landing(clmax_l_designpoint)
        self.tbywto_designpoint = self._tbyw_to_takeoff(clmax_to_designpoint,self.wbysto_designpoint)
        self.wing_area = self.w_to / self.wbysto_designpoint
        self.thrust = self.w_to * self.tbywto_designpoint


if __name__ == "__main__":
    s = Sizing()
    s.display_requirements()
    s.display_setparameter()
    s.calculation()
    s.display_results()
    clmax_l_designpoint = 1.9
    clmax_to_designpoint = 1.7
    s.cal_s_t(clmax_l_designpoint,clmax_to_designpoint)
    s.plot()
