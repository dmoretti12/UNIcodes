#Linear fits

import matplotlib.pyplot as plt
import numpy as np

dofit =  True

k_out2 = [315,332, 345, 348, 389, 433, 433]
# k_out1 = [117, 194, 385,416, 460, 490, 502]
# pkpk1 = [160,163,39,38,41,41,42]
pkpk2 = [53,52,63,48,49,51,46]


icold2 = [4.12, 4.24, 4.31, 4.4, 4.7, 5.02, 5.03]
# icold1 = [2.34, 3.05, 4.73, 5, 5.44, 5.73, 5.73]



pa, = plt.plot(icold2, k_out2,'go--', markersize=3.0)
plt.errorbar(icold2, k_out2, np.sqrt(pkpk2))
plt.xlabel('I_cold (mA)')
plt.ylabel('Koheroen_output (mV)')
plt.title('Finisar pigtailed laser (Rev_3.0, ch2) not potted')

if dofit:
    coef = np.polyfit(icold2, k_out2, 1)
    poly1d_fn = np.poly1d(coef)
    print(coef)
    fittxt='finisar: m= '+str(round(coef[0], 2))+' b= '+str(round(coef[1], 2)) #+ ' RMS= 23.01'
    plt.plot(icold2, poly1d_fn(icold2), label=fittxt)
plt.legend(loc='best')
plt.savefig('finisar pigtailed laser (Rev_3.0, ch2) not potted.png', format = 'png', dpi = 300)

#plt.show()

#RMS = np.sqrt(coef[1]/len(V_noatt))
