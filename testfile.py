from TDMS_Functions import *

def load_tdms_as_Fibers(outpath):
    slx2020=Tdms()
    slx2020.load_variables(outpath)
    slx=Fibres()
    slx.nsamples=mycp(slx2020.nt)
    slx.nchannels=mycp(slx2020.nchannels)
    slx.nt=mycp(slx2020.nt)


    slx.zz=mycp(slx2020.zz)
    slx.data=mycp(slx2020.data)
    slx.tt=mycp(slx2020.tt)
    return slx
