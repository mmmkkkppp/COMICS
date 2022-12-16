import pandas as pd
import astropy.io.fits as fits
import numpy as np
import glob
import os
from scipy.ndimage import gaussian_filter
import pandas as pd
import shutil

## 関数の定義
def readheader(filename):
    return fits.open(filename)[0].header


def readdata(filename):
    return fits.open(filename)[0].data

    return fits.open(filename)[0].data

def q_bsep_posi(input):
    return input[:,::2,:,:]
def q_bsep_nega(input):
    return input[:,1::2,:,:]   
#z方向平均
def mean_z(input):
    return np.mean(input, axis=1, keepdims=True)

#gaussian filter
def gaussfilter(input):
    output = gaussian_filter(input[0, 0, :, :], sigma=6, order=0,
                             output=None, mode='nearest', cval=0.0, truncate=4.0)
    return output


def q_subch(input):
    #z方向については、COMQは1となっているはずである。
    input_dim2 = input[0, 0, :, :]
    #x軸方向に20*240と分割し、z方向に積み上げる

    divstack = np.stack(np.split(input_dim2, 16, axis=1))
    #z方向にmedianをとり、20*240をtileして320*240にする
    beforetile = np.median(divstack, axis=0)  # keepdimsは不要(二次元でいい)
    tilemedian = np.tile(beforetile, 16)
    # print(tilemedian.shape)
    input_to_output = input.copy()
    input_to_output[0, 0, :, :] -= tilemedian
    # input_to_output[0, 0, :, :] = 1
    return input_to_output
#fitsファイル作成用
def makefits(dataname, filename):
    hdu = fits.PrimaryHDU(data=dataname)
    fits.HDUList([hdu]).writeto(str(filename), overwrite=True)

def identify(obsfilename):
    obsheader = readheader(obsfilename)
    #年を取得
    year = obsheader["DATE-OBS"][:4]
    #DARK.csvにアクセス
    darkcsv = r'/mnt/e/'+year +r'/dark/dark.csv'
    darkdf = pd.read_csv(darkcsv)
    #一致条件を指定する
    Q_PIXTIM = obsheader["Q_PIXTIM"]
    Q_RRSTRT = obsheader["Q_RRSTRT"]
    DATE = obsheader["DATE-OBS"]
    Q_YSTRT = obsheader["Q_YSTRT"]
    A = (darkdf["PIXTIM"] == Q_PIXTIM)
    B = (darkdf["RRSTRT"] == Q_RRSTRT)
    C = (darkdf["DATE"] == DATE)
    D = (darkdf["YSTRT"] == Q_YSTRT)
    df = darkdf[A & B & C & D].reset_index(drop=True)
    return df
