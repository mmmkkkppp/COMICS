import subprocess as sp
import sys
import os
from function import *
import pandas as pd
import re
import math

def PSF(file):
    header = readheader(file)["FILTER03"]
    # CとWに挟まれた文字を検索、group(1)でCとWを含めないよう指定、float型に変換
    lambda_ = float(re.search(r'C(.+)W', header).group(1))
    D = 8.08
    PSF_arcsec = 1.028 * lambda_*1e-6 / D * 3600 * 180 / math.pi
    #arcsec からpixelに変換
    PSF_pixel = PSF_arcsec / 0.13
    return PSF_pixel

#sextractorを回す。
#先に、ディレクトリを移動しておく。
filename = sys.argv[1]

data = pd.read_csv('default_.sex', sep='\t', header=None)
#data.iloc[4,0]
#~.catの変更
renamecat = filename.replace('out_obj/','out_obj/cat/').replace('fits', 'cat')#pathをここに入れる
cat = 'CATALOG_NAME     '+ renamecat +'       # name of the output catalog'
data.iloc[4,0] = cat
aftername = 'DETECT_MINAREA   ' + str(PSF(filename))[:5]+'              # min. # of pixels above threshold'
print("DETECT_MINAREA: {}".format(str(PSF(filename))[:5]))
data.iloc[10, 0] = aftername
data.to_csv('finalversion.sex', sep='\t', index=False, header=None)


sp.call(["sex",filename,"-c","finalversion.sex"])
