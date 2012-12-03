from scipy import fftpack, signal
import numpy as np
import math as math
import pylab as py


def FilterDesign((f1, f2), N=100):

    if type(N) is not int:
        raise TypeError, "Filter order must be an integer"

    Ndesign = 256  # assume it is even
    f = fftpack.fftfreq(N, d=0.5)
    af = abs(f)
    desired = np.where((af > f1) & (af < f2), 1, 0)
    hideal = fftpack.ifft(desired).real
    win = signal.get_window('hamming', N)

    # Now we need to get the actual tap coefficients by windowing this result.
    htrunc = fftpack.fftshift(np.r_[hideal[:N / 2], hideal[-(N - 1) / 2:]])
    hfinal = win * htrunc

    return hfinal



def filtfilt(b, a, input_vector):
    """input_vector has shape (n,1)"""
    forward = signal.lfilter(b, a, input_vector, axis=0)
    return np.flipud(signal.lfilter(b, a, np.flipud(forward), axis=0))

def filt(b, a, input_vector):
    """Filter input_vector and correct phase shift"""
    out = signal.lfilter(b, a, np.r_[input_vector, [0] * len(b)])
    return out[len(b) / 2:-len(b) / 2]

def FourierFilter(x, cutoff):
    x_f = fftpack.fft(x)
    n = len(x_f)
    f = np.r_[np.arange((n + 1) / 2), np.arange(n / 2, 0, -1)]
    w = (f < cutoff[1] * n / 2) & (f > cutoff[0] * n / 2)
    return np.real(fftpack.ifft(x_f * w))



def FourierGauss(x, cutoff, sigma, passband_atten=0.5):
    x_f = fftpack.fft(x)
    n = len(x_f)
    try:
        s = sigma[0]
    except TypeError:
        sigma = [sigma, sigma]
    
    #correct cutoff frequencies to take into account transition width
    passband_shift = np.array(sigma) * np.sqrt(-np.log(passband_atten))
    cutoff = [cutoff[0] + passband_shift[0],
             cutoff[1] - passband_shift[1]]
    
    assert cutoff[1] > cutoff[0], ("Left cutoff frequency higher than right. "
                                 "Increase the width of pass-band")
    
    f = np.r_[np.arange((n + 1) / 2), np.arange(n / 2, 0, -1)] * 2. / n
    if sigma[1] > 0:
        w = np.exp(-np.maximum(0, f - cutoff[1]) ** 2 / sigma[1] ** 2)
    else:
        w = 1.*(f < cutoff[1]) 
    if sigma[0] > 0:
        w[f < cutoff[0]] = np.exp(-(f[f < cutoff[0]] - cutoff[0]) ** 2 / sigma[0] ** 2)
    else:
        w[f < cutoff[0]] = 0
    return np.real(fftpack.ifft(x_f * w))

#def FilterResponse(b, a, Fs):
#    """Show frequency and impulse response of a FIR/IIR filter given by
#    coefficients a,b.
#    Note: this function uses phase corrected output vector for calculation of
#    impulse response. See: filt"""
#    w, h = signal.freqz(b, a)
#    
#    N = len(b)
#    t = np.r_[-N / 2:N / 2]
#    y = (t == 0)
#    impz = filt(b, a, y)
#    
#    py.subplot(211); py.semilogy(w / (2 * math.pi) * Fs, abs(h))
#    py.xlabel('Frequency (Hz)')
#    py.subplot(212); py.stem(t, impz); py.xlim((-N / 2, N / 2))
#    py.xlabel('Time (samples)')


   
def remove_base(data):
    """ Removes overall baseline from the data
    function added by: Maja    """
    avg_d = py.mean(data)
    return data - avg_d

def remove_baseloc(data, window):
    """ Removes baseline locally (by defined window size) from the data
    function added by: Maja    """
    weightings = np.repeat(1.0, window) / window
    moved_avg = np.convolve(data, weightings)[window - 1:-(window - 1)]
    
    #import pdb; pdb.set_trace() 
    #print len(moved_avg), len(data)
    m1 = [moved_avg[0]] * (window / 2)
    m2 = [moved_avg[-1]] * (window / (2) - 1)
    moved_avg = np.concatenate((m1, moved_avg, m2))
    if len(moved_avg) < len(data):
        moved_avg = np.hstack([moved_avg, moved_avg[-1]])

    return np.subtract(data, moved_avg), moved_avg

def FilterResponse(b,a,Fs):
    """Show frequency and impulse response of a FIR/IIR filter given by
    coefficients a,b.
    Note: this function uses phase corrected output vector for calculation of
    impulse response. See: filt"""
    w,h=signal.freqz(b,a)

    N=len(b)
    t=plt.r_[-N/2:N/2]
    y=(t==0)
    impz=filt(b,a,y)

    plt.subplot(211); plt.semilogy(w/(2*np.pi)*Fs,abs(h))
    plt.xlabel('Frequency (Hz)')
    plt.subplot(212); plt.stem(t,impz); plt.xlim((-N/2,N/2))
    plt.xlabel('Time (samples)')
    


def bandPass(freq, fs, data, N = 1000):
    # bandpass filter data
  
    b_bp = FilterDesign((freq[0] * 2. / fs, freq[1] * 2. / fs), N) 
    filt_data = filtfilt(b_bp, [1], data)
    #FilterResponse(b_bp,[1],fs)
    return filt_data

def highPass(freq, fs, data, N = 1000):
    # highpass filter data
    b_hp = FilterDesign((freq * 2./ fs, 1), N) 
    fast_data = filtfilt(b_hp, [1], data)   
    return fast_data

def lowPass(freq, fs, data, N = 1000):
    # lowpass filter data 
#    
#    x = np.zeros(1024)
#    x[512]=1
#    Fs = 200 #sampling frequency in Hz
#    N = 100  #filter order 
#    
#    x_fourier = FourierFilter(x, (0.2, 0.6))
#    x_gauss = FourierGauss(x, (0.2, 0.6), (0.1, 0.1) )
   
    
    b_lp = FilterDesign((0, freq * 2. / fs), N) 
    low_data = filtfilt(b_lp, [1], data)
    return low_data

#if __name__ == '__main__':
#
#    x = np.random.randn(1024)
#    Fs = 200 #sampling frequency in Hz
#    N = 100  #filter order 
#    
#    #low pass <20 Hz
#    b_lp = FilterDesign((0, 20. / Fs), N) 
#    x_lp = filtfilt(b_lp, [1], x)
#    
#    #high-pass >50 Hz
#    b_hp = FilterDesign((50. / Fs, 1), N) 
#    x_hp = filtfilt(b_hp, [1], x)
#    
#    #band-pass 20--50 Hz
#    b_bp = FilterDesign((20. / Fs, 50. / Fs), N) 
#    x_bp = filtfilt(b_bp, [1], x)
#
#    import matplotlib.pyplot as plt
#    plt.plot(x_hp, label='highpass')
#    plt.plot(x_lp, label='lowpass')
#    plt.plot(x_bp, label='bandpass')
#    plt.legend()
#    plt.show()
#    
if __name__=='__main__':

    x = np.zeros(1024)
    x[512]=1
    Fs = 200 #sampling frequency in Hz
    N = 100  #filter order 
    
    x_fourier = FourierFilter(x, (0.2, 0.6))
    x_gauss = FourierGauss(x, (0.2, 0.6), (0.1, 0.1) )
    import matplotlib.pyplot as plt
    plt.plot(x_fourier, label='FourierFilter')
    plt.plot(x_gauss, label='FourierGauss')
    plt.legend()
    plt.show()  
    
    
    













