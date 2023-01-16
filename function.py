import pandas as pd
import astropy.io.fits as fits
import numpy as np
import glob
import os
from scipy.ndimage import gaussian_filter
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

#RAWモードの場合にAddモードに変換する
def RAWtoAdd(input, obsheader):
    Q_CHCN = int(obsheader["Q_CHCN"])
    X = np.sum(np.array(np.split(input, Q_CHCN, axis=1)), axis=0)
    obsheader["NAXIS3"] = obsheader["Q_CHCN"]
    obsheader["Q_CHAM"] =0
    return X
    
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
def makefits(dataname, filename, obsfile): 
    hdu = fits.PrimaryHDU(data=dataname,header=fits.getheader(obsfile))
    fits.HDUList([hdu]).writeto(str(filename), overwrite=True)


def obsyear(obs):
    year = readheader(obs)["DATE-OBS"][:4]
    return year

def identify(obsfilename):
    obsheader = readheader(obsfilename)
    # 年を取得
    year = obsheader["DATE-OBS"][:4]
    # DARK.csvにアクセス
    darkcsv = r'/mnt/e/'+year + r'/dark/dark.csv'
    darkdf = pd.read_csv(darkcsv)
    # 一致条件を指定する
    Q_PIXTIM = obsheader["Q_PIXTIM"]
    Q_RRSTRT = obsheader["Q_RRSTRT"]
    DATE = obsheader["DATE-OBS"]
    Q_YSTRT = obsheader["Q_YSTRT"]
    NAXIS4 = obsheader["NAXIS4"]
    A = (darkdf["PIXTIM"] == Q_PIXTIM)
    B = (darkdf["RRSTRT"] == Q_RRSTRT)
    C = (darkdf["DATE"] == DATE)
    D = (darkdf["YSTRT"] == Q_YSTRT)
    E = (darkdf["NAXIS4"] == NAXIS4)
    df = darkdf[A & B & C & D & E].reset_index(drop=True)
    return df



def make1Expdark(obsfile):
    df = identify(obsfile)
    year = obsyear(obsfile)
    darkbasefolder = r'/mnt/e/' + year + r'/dark/'
    obsdata = readdata(obsfile)
    # obsdataと同じshapeの0のnumpyを作成
    darkdata = np.zeros(list(obsdata.shape))
    if len(df) == 0:
        # print("darkが存在しない")
        return None
        # ここにはその場合の処理を記載する。

    # 特定したダークファイル一枚一枚に対して、Q_CHEBが同じものを抽出??

    # df[df["Q_CHEB"].duplicated()]
    for darkfile in df["FRAME_ID"]:
        darkfilepath = darkbasefolder + darkfile + r'.fits'
        darkdata += np.mean(readdata(darkfilepath), axis=1, keepdims=True)
    # meanをとる
    meandark = darkdata / len(df)
    # 1Expあたりに直す
    meandark_1 = meandark / int(df["Q_CHEB"][0])
    return meandark_1

#積分時間をとってくる関数
def exptime(obsfile):
    obsheader = readheader(obsfile)
    Q_1FRAME = obsheader["Q_1FRAME"]
    Q_CHCN = obsheader["Q_CHCN"]
    print(Q_1FRAME)
    print(Q_CHCN)
    exptime = Q_1FRAME * Q_CHCN / 2
    return exptime

def exptime1(obsfile):
    obsheader = readheader(obsfile)
    Q_1FRAME = obsheader["Q_1FRAME"]
    Q_CHCN = obsheader["Q_CHCN"]
    print(Q_1FRAME)
    print(Q_CHCN)
    exptime1 = int(Q_1FRAME) * int(Q_CHCN) / 2
    return exptime1


def make_obj(obsfile):
    
# dark_1を取得
    
    dark_1 = make1Expdark(obsfile)
    if dark_1 is None:
        return None
    # Q_CHEB倍する
    obsheader = readheader(obsfile)
    Q_CHEB = obsheader["Q_CHEB"]
    dark_CHEB = dark_1 * int(Q_CHEB)
    # obsdata
    # この時点で、RAWをAddに変換する。
    if (obsheader["Q_CHAM"] ==0):
        obsdata_RAW = readdata(obsfile)
        Q_CHCN = obsheader["Q_CHCN"]
        A = np.array(np.split(obsdata_RAW, Q_CHCN, axis=1))
        obsdata = np.sum(A, axis=0)
    elif (obsheader["Q_CHAM"]==1):
        obsdata = readdata(obsfile)
    elif (obsheader["Q_CHAM"]==2):
        print("ECOモードによりskip")
    
    
    # skydata作成
    skydata = obsdata - dark_CHEB
    # q_bsepによりposiとnegaにわける
    # Addモードの場合だけでやってみる。
    skydata_p = q_bsep_posi(skydata)
    skydata_n = q_bsep_nega(skydata)

    # z方向平均
    sky_pa = mean_z(skydata_p)
    sky_na = mean_z(skydata_n)

    # gaussian　平均
    sky_paG = gaussfilter(sky_pa)
    sky_naG = gaussfilter(sky_na)

    # Flatの作成
    sky_paF = sky_pa / sky_paG
    sky_naF = sky_na / sky_naG
    
    comqfile = obsfile.replace("COMA", "COMQ")

    comq_obs = readdata(comqfile)
    obj_obs = q_subch(comq_obs)
    # Flatで割る
    obj_obs_datP0 = obj_obs / sky_naF
    obj_obs_datN0 = obj_obs / sky_paF
    
    #Exptimeで割る
    obj_obs_datP0_per1sec = obj_obs_datP0 / exptime(comqfile)
    obj_obs_datN0_per1sec = obj_obs_datN0 / exptime(comqfile)
    
    os.makedirs(r'/mnt/e/'+obsyear(obsfile)+r'/out_obj/',exist_ok=True)
    outputfol_P = r'/mnt/e/'+obsyear(obsfile)+r'/out_obj/obj_'+obsfile[-10:-5]+r'_datP0_per1sec.fits'
    makefits(obj_obs_datP0_per1sec,outputfol_P, comqfile)
    
    outputfol_N = r'/mnt/e/' + \
        obsyear(obsfile)+r'/out_obj/obj_' + \
        obsfile[-10:-5]+r'_datN0_per1sec.fits'
    makefits(obj_obs_datN0_per1sec, outputfol_N, comqfile)
    
    return obj_obs_datP0_per1sec, obj_obs_datN0_per1sec
