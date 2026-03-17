import astropy
from astropy.table import Table, join, vstack
from astropy.time import Time
from astropy.io import fits
from astroquery.jplhorizons import Horizons
import mastcasjobs
from PIL import Image, UnidentifiedImageError
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
import time
import requests
import os
import tkinter as tk
from tkinter import simpledialog, messagebox


# 配置选项
saveplots = True  # 是否保存生成的图表
save_data = True  # 是否将观测表和图像保存到本地
use_jpl = True  # 是否尝试从 JPL Horizons 获取高精度星历（设为False跳过）

astropy.conf.max_width = 150


if not os.environ.get('CASJOBS_USERID') or not os.environ.get('CASJOBS_PW'):
    root = tk.Tk()
    root.withdraw()  # Hide main window

    userid = simpledialog.askstring("CasJobs Login", "Enter Casjobs UserID:")
    if userid is None or userid.strip() == "":
        messagebox.showerror("Login Failed", "Username is required.")
        raise SystemExit("User cancelled or did not provide username.")

    password = simpledialog.askstring("CasJobs Login", "Enter Casjobs password:", show='*')
    if password is None or password.strip() == "":
        messagebox.showerror("Login Failed", "Password is required.")
        raise SystemExit("User cancelled or did not provide password.")

    os.environ['CASJOBS_USERID'] = userid.strip()
    os.environ['CASJOBS_PW'] = password.strip()

    root.destroy()

try:
    jobs = mastcasjobs.MastCasJobs(context="MyDB")
    print("Confirmed MAST CasJobs account access")
except Exception as e:
    print(f"Warning: Failed to connect to CasJobs: {e}")
    print("Continuing without CasJobs (may be OK if not needed).")



def cadc_ssos_query(object_name, search="bynameCADC", epoch1=54985, epoch2=57079,
                    xyres="no", telinst="Pan-STARRS1", lang="en", format="tsv",
                    url="https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl"):
    """使用CADC移动天体查询功能按目标名称查找PS1观测数据。历元参数给出了PS1观测数据的MJD范围。
    唯一可能有修改价值的参数是"search"。例如，使用search="bynameHorizons"，可以借助JPL Horizons星历表进行搜索，而非使用CADC缓存的星历表。
    若CADC的缓存已过期，这种方式可能会派上用场。这将返回一个包含观测数据的astropy表。
    """
    r = requests.get(url, params=dict(lang=lang, object=object_name, search=search,
                                      epoch1=epoch1, epoch2=epoch2, xyres=xyres,
                                      telinst=telinst, format=format))
    r.raise_for_status()
    # Check for "not found"
    msg = r.text.split('\n', 1)[0]
    if "not found" in msg.lower():
        if msg.endswith("<br/>"):
            msg = msg[:-5]
        raise ValueError(f"Error from CADC search: {msg}")
    tab = Table.read(r.text, format="ascii.csv", delimiter="\t",
                     converters=dict(Image_target=str))
    return tab



def getimages(tab, format="jpg", verbose=False, **kw):
    """返回由包含常规列的表格指定的图像列表
    （常规列包括：projcell、subcell、image、赤经(ra)和赤纬(dec)列）
    参数说明：
    format = 数据格式（可选值为 "jpg"、"png"、"fits"）
    verbose = 若为真，则打印进度信息
    传递给 tab2cutout 函数的其他参数：
    size = 提取图像的尺寸（单位：像素，像素分辨率为 0.25 角秒/像素）
    output_size = 输出（展示）图像的尺寸（单位：像素，默认值 = size）。
                  对于 fits 格式的图像，output_size 参数无作用。
    racol、deccol、imagecol = 表格列的名称
    返回值：以 PIL Image 图像对象 或 astropy fits 库的 hdulist 列表对象形式返回图像
    """
    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be jpg, png or fits")
    urls = tab2cutout(tab, format=format, verbose=verbose, **kw)
    images = []
    if verbose:
        t0 = time.time()
        print(f"{time.time() - t0:.1f} s: completed 0 of {len(urls)} images", end="")
    for i, url in enumerate(urls, start=1):
        r = requests.get(url)
        if format == "fits":
            images.append(fits.open(BytesIO(r.content)))
        else:
            try:
                images.append(Image.open(BytesIO(r.content)))
            except UnidentifiedImageError as e:
                print(f"\nImage not found {url}")
                images.append(None)
        if verbose:
            print(f"\r{time.time() - t0:.1f} s: completed {i} of {len(urls)} images", end="")
    if verbose:
        print()
    return images


def tab2path(tab, imagecol='Image', verbose=False):
    """将包含常规列（projcell、subcell、filter、expStart）的表格转换为 warp滤镜路径，返回文件路径列表"""
    for c in ('projcell', 'subcell', imagecol):
        assert c in tab.colnames, f"'{c}' column not found in table"
    rv = []
    for i, row in enumerate(tab):
        projcell = f"{row['projcell']:04d}"
        subcell = f"{row['subcell']:03d}"
        image = row[imagecol]
        rv.append(f"/rings.v3.skycell/{projcell}/{subcell}/{image}.fits")
    return rv


def tab2cutout(tab, size=64, output_size=256, format="jpg", autoscale=99.5, asinh=True,
               racol="ra", deccol="dec", imagecol="Image", verbose=False,
               urlbase="https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"):
    """将包含常规列的表格转换为 warp 图像裁剪对应的 URL 地址"""
    for c in ('projcell', 'subcell', imagecol, racol, deccol):
        assert c in tab.colnames, f"'{c}' column not found in table"
    paths = tab2path(tab, imagecol=imagecol, verbose=verbose)
    rv = []
    for row, path in zip(tab, paths):
        params = dict(red=path, ra=row[racol], dec=row[deccol], size=size, format=format)
        if format != "fits":
            params['autoscale'] = autoscale
            if output_size is not None:
                params['output_size'] = output_size
            if asinh:
                params['asinh'] = 1
        sparams = '&'.join([f"{key}={str(value)}" for key, value in params.items()])
        rv.append(f"{urlbase}?{sparams}")
    return rv


# pixel scale is 0.25 arcsec
pixscale = 0.25

# table of rings info
rings = Table.read("""zone projcell nband dec dec_min dec_max xcell ycell crpix1 crpix2
13 487 72 -38.0 -39.98569922229713 -35.98748871852416 6305 6279 239.5 238.0
14 559 76 -34.0 -35.98748871852416 -31.988947843740302 6265 6274 237.5 237.5
15 635 79 -30.0 -31.988947843740302 -27.990585819085737 6274 6269 239.5 241.5
16 714 82 -26.0 -27.990585819085737 -23.991856317278298 6255 6265 241.5 240.5
17 796 84 -22.0 -23.991856317278298 -19.99333106003115 6279 6261 241.0 242.0
18 880 86 -18.0 -19.99333106003115 -15.994464100667857 6274 6258 241.5 238.0
19 966 88 -14.0 -15.994464100667857 -11.995671545879146 6242 6254 240.5 239.5
20 1054 89 -10.0 -11.995671545879146 -7.997014822849273 6248 6250 240.0 242.0
21 1143 89 -6.0 -7.997014822849273 -3.998086931795508 6291 6247 237.5 240.0
22 1232 90 -2.0 -3.998086931795508 1.3877787807814457e-17 6240 6243 240.0 242.0
23 1322 90 2.0 1.3877787807814457e-17 3.998086931795508 6240 6243 240.0 242.0
24 1412 89 6.0 3.998086931795508 7.997014822849273 6291 6247 237.5 240.0
25 1501 89 10.0 7.997014822849273 11.995671545879146 6248 6250 240.0 242.0
26 1590 88 14.0 11.995671545879146 15.994464100667857 6242 6254 240.5 239.5
27 1678 86 18.0 15.994464100667857 19.99333106003115 6274 6258 241.5 238.0
28 1764 84 22.0 19.99333106003115 23.991856317278298 6279 6261 241.0 242.0
29 1848 82 26.0 23.991856317278298 27.990585819085737 6255 6265 241.5 240.5
30 1930 79 30.0 27.990585819085737 31.988947843740302 6274 6269 239.5 241.5
31 2009 76 34.0 31.988947843740302 35.98748871852416 6265 6274 237.5 237.5
32 2085 72 38.0 35.98748871852416 39.98569922229713 6305 6279 239.5 238.0
33 2157 68 42.0 39.98569922229713 43.98384988850235 6320 6284 239.5 240.0
34 2225 64 46.0 43.98384988850235 47.98179099058205 6307 6289 238.0 242.0
35 2289 60 50.0 47.98179099058205 51.97972846465295 6261 6295 241.0 239.5
36 2349 55 54.0 51.97972846465295 55.97676677254534 6283 6302 239.0 242.0
37 2404 50 58.0 55.97676677254534 59.97411834356914 6278 6311 238.5 237.5
38 2454 45 62.0 59.97411834356914 63.97005521318337 6240 6319 240.0 241.0
39 2499 39 66.0 63.97005521318337 67.96475284198297 6307 6333 239.5 240.0
40 2538 33 70.0 67.96475284198297 71.95817764373305 6365 6350 238.5 240.0
41 2571 27 74.0 71.95817764373305 75.94974140116577 6413 6371 240.5 242.0
42 2598 21 78.0 75.94974140116577 79.93968190571101 6452 6398 240.0 242.0
43 2619 15 82.0 79.93968190571101 83.93610604798963 6481 6429 241.0 242.0
44 2634 9 86.0 83.93610604798963 87.96940975734594 6501 6419 239.0 239.5
45 2643 1 90.0 87.96940975734594 90.0001 6240 6240 239.5 240.0
""", format="ascii.csv", delimiter=" ")
dec_limit = rings['dec_min'].min()



def findskycell(ra, dec):
    if np.isscalar(ra) and np.isscalar(dec):
        return _findskycell_array(np.array([ra]), np.array([dec]))
    if len(ra) == len(dec):
        return _findskycell_array(np.asarray(ra), np.asarray(dec))
    else:
        raise ValueError("ra and dec must both be scalars or be matching length arrays")


def _findskycell_array(ra, dec):
    if ra.ndim != 1 or dec.ndim != 1:
        raise ValueError("ra and dec must be 1-D arrays")
    index = np.arange(len(ra), dtype=int)
    idec = np.searchsorted(rings['dec_max'], dec)
    nearpole = np.where(idec >= len(rings) - 2)
    idec[nearpole] = len(rings) - 2
    nband = rings['nband'][idec]
    nra = ra % 360.0
    ira = np.rint(nra * nband / 360.0).astype(int) % nband
    projcell = rings['projcell'][idec] + ira
    dec_cen = rings['dec'][idec]
    ra_cen = ira * 360.0 / nband
    x, y = sky2xy_tan(nra, dec, ra_cen, dec_cen)
    pad = 480
    if len(nearpole[0]) > 0:
        projcell2 = rings['projcell'][-1]
        dec_cen2 = rings['dec'][-1]
        ra_cen2 = 0.0
        x2, y2 = sky2xy_tan(nra[nearpole], dec[nearpole], ra_cen2, dec_cen2)
        use2 = poleselect(x[nearpole], y[nearpole], x2, y2, rings[-2], rings[-1], pad)
        if use2.any():
            wuse2 = np.where(use2)[0]
            w2 = tuple(xx[wuse2] for xx in nearpole)
            idec[w2] = len(rings) - 1
            nband[w2] = 1
            ira[w2] = 0
            projcell[w2] = projcell2
            dec_cen[w2] = dec_cen2
            ra_cen[w2] = ra_cen2
            x[w2] = x2[wuse2]
            y[w2] = y2[wuse2]
    px = rings['xcell'][idec] - pad
    py = rings['ycell'][idec] - pad
    k = np.rint(4.5 + x / px).astype(int).clip(0, 9)
    j = np.rint(4.5 + y / py).astype(int).clip(0, 9)
    subcell = 10 * j + k
    crpix1 = rings['crpix1'][idec] + px * (5 - k)
    crpix2 = rings['crpix2'][idec] + py * (5 - j)
    ximage = x + crpix1
    yimage = y + crpix2
    iring = idec
    w = np.where(dec < dec_limit)
    projcell[w] = 0
    subcell[w] = 0
    crpix1[w] = 0
    crpix2[w] = 0
    ximage[w] = 0
    yimage[w] = 0
    return Table([ra, dec, index, projcell, subcell, ra_cen, dec_cen, crpix1, crpix2, ximage, yimage, iring],
                 names="ra,dec,index,projcell,subcell,crval1,crval2,crpix1,crpix2,x,y,iring".split(","))


def poleselect(x1, y1, x2, y2, rings1, rings2, pad):
    nx1 = 10 * (rings1['xcell'] - pad) + pad
    ny1 = 10 * (rings1['ycell'] - pad) + pad
    nx2 = 10 * (rings2['xcell'] - pad) + pad
    ny2 = 10 * (rings2['ycell'] - pad) + pad
    d1 = np.minimum(np.minimum(x1 + nx1 // 2, nx1 // 2 - 1 - x1), np.minimum(y1 + ny1 // 2, ny1 // 2 - 1 - y1))
    d2 = np.minimum(np.minimum(x2 + nx2 // 2, nx2 // 2 - 1 - x2), np.minimum(y2 + ny2 // 2, ny2 // 2 - 1 - y2))
    return d1 < d2


def sky2xy_tan(ra, dec, ra_cen, dec_cen, crpix=(0.0, 0.0)):
    dtor = np.pi / 180
    cd00 = -pixscale * dtor / 3600
    cd01 = 0.0
    cd10 = 0.0
    cd11 = -cd00
    determ = cd00 * cd11 - cd01 * cd10
    cdinv00 = cd11 / determ
    cdinv01 = -cd01 / determ
    cdinv10 = -cd10 / determ
    cdinv11 = cd00 / determ
    cos_crval1 = np.cos(dtor * dec_cen)
    sin_crval1 = np.sin(dtor * dec_cen)
    radif = (ra - ra_cen) * dtor
    w = np.where(radif > np.pi)
    radif[w] -= 2 * np.pi
    w = np.where(radif < -np.pi)
    radif[w] += 2 * np.pi
    decrad = dec * dtor
    cos_dec = np.cos(decrad)
    sin_dec = np.sin(decrad)
    cos_radif = np.cos(radif)
    sin_radif = np.sin(radif)
    h = sin_dec * sin_crval1 + cos_dec * cos_crval1 * cos_radif
    xsi = cos_dec * sin_radif / h
    eta = (sin_dec * cos_crval1 - cos_dec * sin_crval1 * cos_radif) / h
    xdif = cdinv00 * xsi + cdinv01 * eta
    ydif = cdinv10 * xsi + cdinv11 * eta
    return (xdif + crpix[0], ydif + crpix[1])


def xy2sky_tan(x, y, ra_cen, dec_cen, crpix=(0.0, 0.0)):
    dtor = np.pi / 180
    cd00 = -pixscale * dtor / 3600
    cd01 = 0.0
    cd10 = 0.0
    cd11 = -cd00
    cos_crval1 = np.cos(dtor * dec_cen)
    sin_crval1 = np.sin(dtor * dec_cen)
    xdif = x - crpix[0]
    ydif = y - crpix[1]
    xsi = cd00 * xdif + cd01 * ydif
    eta = cd10 * xdif + cd11 * ydif
    beta = cos_crval1 - eta * sin_crval1
    ra = np.arctan2(xsi, beta) + dtor * ra_cen
    gamma = np.sqrt(xsi ** 2 + beta ** 2)
    dec = np.arctan2(eta * cos_crval1 + sin_crval1, gamma)
    return (ra / dtor, dec / dtor)



# 下载原始 warp FITS 文件
def download_warp_fits(tab, download_dir, base_url="https://ps1images.stsci.edu/data/ps1/data/"):
    """
    下载表格中每一行对应的原始 warp FITS 文件。
    使用 tab2path 构建路径，并从 PS1 数据服务器下载。
    """
    os.makedirs(download_dir, exist_ok=True)
    paths = tab2path(tab, imagecol='Image')  # 获取相对路径
    downloaded_files = []

    for i, (row, path) in enumerate(zip(tab, paths)):
        filename = os.path.basename(path) + ".fits"
        url = base_url + path.lstrip("/") + ".fits"
        filepath = os.path.join(download_dir, filename)

        # 跳过已存在的文件
        if os.path.exists(filepath):
            print(f"  已存在，跳过: {filepath}")
            downloaded_files.append(filepath)
            continue

        print(f"[{i + 1}/{len(tab)}] 下载: {url}")
        try:
            r = requests.get(url, stream=True, timeout=30)
            r.raise_for_status()
            with open(filepath, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f" 成功保存: {filepath}")
            downloaded_files.append(filepath)
        except requests.exceptions.RequestException as e:
            print(f" 下载失败 {url}: {e}")
        except Exception as e:
            print(f" 保存失败 {filepath}: {e}")

    return downloaded_files



source = "2016 OR17"  # 目标编号
output_dir_base = f"PS1_data_{source}"
os.makedirs(output_dir_base, exist_ok=True)

# 查询 CADC
cadc_tab = cadc_ssos_query(source)
print(f"Returned table has {len(cadc_tab)} rows")
print(cadc_tab[:8])

# 初始化位置列
ra_col = cadc_tab['Object_RA'].copy()
dec_col = cadc_tab['Object_Dec'].copy()
fullname = f"Object {source}"

# 尝试从 JPL 获取高精度位置
if use_jpl:
    try:
        print("Attempting to retrieve high-precision ephemerides from JPL Horizons...")
        jdepochs = Time(cadc_tab["MJD"], scale="tai", format="mjd").utc.jd.tolist()
        tlist = []
        block = 50
        for k in range(0, len(jdepochs), block):
            obj = Horizons(id=source, location='F51', epochs=jdepochs[k:k + block])
            tlist.append(obj.ephemerides(extra_precision=True))
        jpltab = vstack(tlist)
        fullname = jpltab['targetname'][0]
        print(f"Got positions from Horizons for {len(jdepochs)} epochs with {len(jpltab)} rows for {fullname}")
        assert np.allclose(jdepochs, jpltab['datetime_jd'].data)
        t2 = Time(jpltab['datetime_jd'].data, format='jd', scale='utc').tai.mjd
        assert np.allclose(t2, cadc_tab["MJD"].data)
        ra_col = jpltab['RA']
        dec_col = jpltab['DEC']
        cadc_tab['illumination'] = jpltab['illumination']
        cadc_tab['illumination'].format = ".2f"
        cadc_tab['delta'] = jpltab['delta']
        cadc_tab['delta'].format = ".2f"
        cadc_tab['delta'].unit = "AU"
    except Exception as e:
        print(f"Warning: Could not retrieve data from JPL Horizons: {e}")
        print("Falling back to CADC-provided positions.")

# 更新主表中的位置
for c in ('ra', 'dec'):
    if c in cadc_tab.colnames:
        cadc_tab[c] = [ra_col, dec_col][['ra', 'dec'].index(c)]
    else:
        cadc_tab.add_column([ra_col, dec_col][['ra', 'dec'].index(c)], name=c,
                            index=cadc_tab.colnames.index('Object_RA'))

# 计算偏差
cadc_tab['dra'] = (cadc_tab['ra'] - cadc_tab['Object_RA']) * 3600 * np.cos(np.radians(cadc_tab['dec']))
cadc_tab['ddec'] = (cadc_tab['dec'] - cadc_tab['Object_Dec']) * 3600
for c in ('ra', 'dec', 'Object_RA', 'Object_Dec'):
    cadc_tab[c].format = ".6f"
for c in ('dra', 'ddec'):
    cadc_tab[c].format = ".2f"
    cadc_tab[c].unit = "arcsec"

print(
    f"dra range  {cadc_tab['dra'].min():.2f} {cadc_tab['dra'].max():.2f} MAD {np.median(np.abs(cadc_tab['dra'])):.2f} arcsec")
print(
    f"ddec range {cadc_tab['ddec'].min():.2f} {cadc_tab['ddec'].max():.2f} MAD {np.median(np.abs(cadc_tab['ddec'])):.2f} arcsec")

# 添加 sky cell 信息
rv = findskycell(cadc_tab['ra'], cadc_tab['dec'])
for c in rv.colnames:
    if c not in cadc_tab.colnames:
        cadc_tab[c] = rv[c]
cadc_tab["crval1"].format = ".6f"
for c in ("x", "y"):
    cadc_tab[c].format = ".2f"

# 处理重复项
umjd, umjd_counts = np.unique(cadc_tab['MJD'], return_counts=True)
nuniq = len(umjd)
print(f"cadc_tab has {len(cadc_tab) - nuniq} duplicate entries, leaving {nuniq} unique entries")

keep = np.array([image.startswith(f"rings.v3.skycell.{projcell:04d}.{subcell:03d}.")
                 for image, projcell, subcell in
                 zip(cadc_tab['Image'].data, cadc_tab['projcell'].data, cadc_tab['subcell'].data)])

if keep.sum() != nuniq:
    jtest = join(cadc_tab['MJD', 'Image', 'projcell', 'subcell'], cadc_tab['MJD', 'Image_target'][keep],
                 join_type='left')
    w = np.where(jtest['Image_target'].mask)[0]
    if len(w) > nuniq - keep.sum():
        j = len("rings.v3.skycell.nnnn.nnn.")
        cmulti = {}
        nchanged = 0
        for i in w:
            tt = jtest['MJD'][i]
            cmulti[tt] = cmulti.get(tt, 0) + 1
            if cmulti[tt] == 2:
                projcell = jtest['projcell'][i]
                subcell = jtest['subcell'][i]
                image = jtest['Image'][i]
                kk = np.where((cadc_tab['MJD'] == tt) &
                              (cadc_tab['projcell'] == projcell) &
                              (cadc_tab['subcell'] == subcell) &
                              (cadc_tab['Image'] == image))[0]
                assert len(kk) == 1
                kk = kk[0]
                oldname = cadc_tab['Image'][kk]
                if oldname.startswith(f"rings.v3.skycell.{projcell:04d}."):
                    cadc_tab['Image'][kk] = f"rings.v3.skycell.{projcell:04d}.{subcell:03d}.{oldname[j:]}"
                    nchanged += 1
        if nchanged > 0:
            keep = np.array([image.startswith(f"rings.v3.skycell.{projcell:04d}.{subcell:03d}.")
                             for image, projcell, subcell in
                             zip(cadc_tab['Image'].data, cadc_tab['projcell'].data, cadc_tab['subcell'].data)])
            jtest = join(cadc_tab['MJD', 'Image', 'projcell', 'subcell'], cadc_tab['MJD', 'Image_target'][keep],
                         join_type='left')
            w = np.where(jtest['Image_target'].mask)[0]
    assert len(w) == nuniq - keep.sum()
    assert (cadc_tab['Image'][w] == jtest['Image'][w]).all()
    new_projcell = [int(x.split('.')[0]) for x in cadc_tab['Image_target'][w].data]
    new_subcell = [int(x.split('.')[1]) for x in cadc_tab['Image_target'][w].data]
    print(f"Patching projcell and subcell for {len(w)} images")
    cadc_tab['projcell'][w] = new_projcell
    cadc_tab['subcell'][w] = new_subcell
    keep = np.array([image.startswith(f"rings.v3.skycell.{projcell:04d}.{subcell:03d}.")
                     for image, projcell, subcell in
                     zip(cadc_tab['Image'].data, cadc_tab['projcell'].data, cadc_tab['subcell'].data)])
    assert keep.sum() == nuniq
    print("Modified lines:")
    print(cadc_tab['Image', 'MJD', 'Filter', 'Image_target', 'projcell', 'subcell'][w])

cadc_tab = cadc_tab[keep]

# 时间范围
emin = Time(54985.0, format='mjd')
emax = Time(57079.0, format='mjd')
print(f"Searching for {source} from {emin.isot} to {emax.isot} ({(emax.mjd - emin.mjd):.0f} days)")

# 绘图：轨道轨迹
plt.rcParams.update({"font.size": 14})
plt.figure(1, (12, 6))

ratab = cadc_tab['ra']
dec_plot = cadc_tab['dec']
jpl_mjd = cadc_tab['MJD']

if use_jpl:
    try:
        obj = Horizons(id=source, location='F51', epochs=dict(start=emin.isot, stop=emax.isot, step='1d'))
        full_jpltab = obj.ephemerides(extra_precision=True)
        jpl_mjd = Time(full_jpltab['datetime_jd'], format='jd').mjd
        ratab = full_jpltab['RA'].data.copy()
        dec_plot = full_jpltab['DEC'].data.copy()
        # Unwrap RA
        w = np.where(ratab[1:] > ratab[:-1] + 180)[0]
        for k in w: ratab[k + 1:] -= 360
        w = np.where(ratab[1:] < ratab[:-1] - 180)[0]
        for k in w: ratab[k + 1:] += 360
        while ratab.min() < 0: ratab += 360
        plt.plot(ratab, dec_plot, '-', alpha=0.6, color='gray', label="Orbital path")
        mstep = max(1, int(30 / (jpl_mjd[1] - jpl_mjd[0]) + 0.5))
        plt.plot(ratab[::mstep], dec_plot[::mstep], 'o', color="tab:blue", label="Monthly positions")
    except Exception as e:
        print(f"Could not get full ephemeris: {e}")

ra_interp = np.interp(cadc_tab['MJD'], jpl_mjd, ratab)
dec_interp = np.interp(cadc_tab['MJD'], jpl_mjd, dec_plot)
plt.plot(ra_interp, dec_interp, 'o', color="red", label=f"{len(cadc_tab)} PS1 obs")
plt.xlabel("Unwrapped RA [deg]")
plt.ylabel("Dec [deg]")
plt.legend(loc="best")
plt.title(f"{fullname} from {emin.isot[:10]} to {emax.isot[:10]}")
plt.tight_layout()
if saveplots:
    plt.savefig(os.path.join(output_dir_base, "sky_path.png"), facecolor="white", dpi=150)

# 下载图像（缩略图）
imsize = 65
output_size = 130
images = getimages(cadc_tab, size=imsize, output_size=output_size, verbose=True)
print(f"Got {len(images)} images")

# 过滤空白图像
goodimages = [img for img in images if img is not None and np.array(img).min() != 255]
igood = [i for i, img in enumerate(images) if img is None or np.array(img).min() != 255]
ngood = len(goodimages)
goodtab = cadc_tab[igood]
goodtab.add_column(np.arange(len(goodtab)), name="SrcID", index=0)
print(f"{ngood} of {len(images)} are not blank")


# 查询PS1目录数据库获取检测结果
dtab = None
try:
    print("\nQuerying PS1 catalog for detections...")
    cnames = ["SrcID", "ra", "dec", "Filter", "MJD"]
    with StringIO() as fh:
        goodtab[cnames].write(fh, format="ascii.csv")
        hdata = fh.getvalue()
    table = "jplhorizons_temp"
    jobs = mastcasjobs.MastCasJobs(context="MyDB", request_type="POST")
    jobs.drop_table_if_exists(table)
    print(f"Dropped MyDB.{table}")
    jobs.upload_table(table, hdata, exists=False)
    print(f"Uploaded table to MyDB.{table}")

    query = f"""
    select t.SrcID, dbo.fDistanceArcminEq(d.ra,d.dec,t.ra,t.dec)*60 as darcsec,
           t.ra as hra, t.dec as hdec,
           t.filter as hfilter, t.mjd as hmjd,
           d.objID, d.detectID, d.filterID, d.imageID, d.obsTime, d.ra, d.dec,
           d.raErr, d.decErr, d.expTime, d.psfFlux, d.psfFluxErr,
           d.psfMajorFWHM, d.psfMinorFWHM, d.psfTheta, d.psfCore,
           d.psfQfPerfect, d.psfChiSq, d.apFlux, d.apFluxErr,
           d.infoFlag, d.infoFlag2, d.infoFlag3
    from MyDB.{table} t
    cross apply fGetNearbyObjEq(t.ra, t.dec, 2.0/60.0) nb
    join Detection d on nb.objid=d.objid 
        and d.obsTime between t.mjd-0.5*d.expTime/86400 and t.mjd+0.5*d.expTime/86400
    order by d.obsTime, darcsec
    """

    dtab = jobs.quick(query, context="PanSTARRS_DR2")
    print(f"Retrieved {len(dtab)} row table")
    dtab['filterID'] = dtab['filterID'].astype(np.uint8)
    dtab['darcsec'].format = ".2f"

    jobs.drop_table(table)

    # 显示位置差异的散点图
    plt.figure(3, (12, 12))
    wgood = np.where(dtab['psfQfPerfect'] > 0.95)[0]
    wbad = np.where(dtab['psfQfPerfect'] <= 0.95)[0]
    dra = (dtab['ra'] - dtab['hra']) * 3600 * np.cos(np.radians(dtab['hdec']))
    ddec = (dtab['dec'] - dtab['hdec']) * 3600
    plt.plot(dra[wgood], ddec[wgood], 'o', color='lightblue', label=f"{len(wgood)} good points")
    plt.plot(dra[wbad], ddec[wbad], 'o', color='red', label=f"{len(wbad)} bad points")
    plt.xlabel("Delta RA [arcsec]")
    plt.ylabel("Delta Dec [arcsec]")
    plt.gca().set_aspect(1.0)
    plt.legend(loc="best", title="psfQfPerfect flags")
    plt.title(f"PS1 catalog detections for {fullname}")
    if saveplots:
        plt.savefig(os.path.join(output_dir_base, "position-diffs.png"), facecolor="white", bbox_inches="tight")

except Exception as e:
    print(f"Warning: Could not query PS1 catalog: {e}")
    print("Continuing without source position markings.")


# 显示并保存图像网格（新增源位置标记）
# 显示并保存图像网格，标记源位置
ncols = 7
nrows = (ngood + ncols - 1) // ncols
hsize = ncols * 2.5
vsize = nrows * 2.5
extent = np.array([1, -1, -1, 1]) * 0.25 * imsize / 2

plt.rcParams.update({"font.size": 10})
plt.figure(2, (hsize, vsize))

# 如果有检测结果，准备标记数据
goodmatch = None
if dtab is not None and len(dtab) > 0:
    # 将 goodtab 表与检测表关联
    goodmatch = join(goodtab, dtab['SrcID', 'ra', 'dec', 'psfFlux', 'psfQfPerfect'], keys='SrcID', join_type='left')

    # 计算源目标相对于图像中心的偏移量
    dx = (goodmatch['ra_2'] - goodmatch['ra_1']) * 3600 * np.cos(np.radians(goodmatch['dec_1']))
    dy = (goodmatch['dec_2'] - goodmatch['dec_1']) * 3600

    # 将流量值转换为 AB 星等
    mask = ~np.isnan(goodmatch['psfFlux'])
    mag = np.full(len(goodmatch), np.nan)
    mag[mask] = -2.5 * np.log10(goodmatch['psfFlux'][mask]) + 8.90
    goodmatch['mag'] = mag

for i in range(ngood):
    plt.subplot(nrows, ncols, i + 1)
    plt.imshow(goodimages[i], origin="upper", cmap="gray", vmin=0, vmax=255, extent=extent)

    # 若源目标位置可用，则标记该位置
    if goodmatch is not None and not np.isnan(goodmatch['ra_2'][i]) and not np.isnan(goodmatch['dec_2'][i]):
        # Calculate offset from image center
        dx_val = (goodmatch['ra_2'][i] - goodmatch['ra_1'][i]) * 3600 * np.cos(np.radians(goodmatch['dec_1'][i]))
        dy_val = (goodmatch['dec_2'][i] - goodmatch['dec_1'][i]) * 3600

        # 根据质量确定颜色
        if goodmatch['psfQfPerfect'][i] > 0.95:
            color = "lightblue"
        else:
            color = "red"

        # 用对应颜色绘制源目标的位置
        plt.plot(dx_val, dy_val, 'o', color=color, markerfacecolor="none", markersize=15, markeredgewidth=2)

        if not np.isnan(goodmatch['mag'][i]):
            plt.text(dx_val, dy_val - 5, f"{goodmatch['mag'][i]:.2f}", color="yellow", fontsize=8,
                     ha='center', va='top', weight='bold')

    plt.title(f"{goodtab['Filter'][i]} {goodtab['MJD'][i]:.5f}")
plt.tight_layout()
plt.suptitle(f"{ngood} {imsize}x{imsize} pixel ({imsize / 4:.0f}x{imsize / 4:.0f} arcsec) PS1 images for {fullname}",
             y=1.01, fontsize=16, fontweight="bold")
if saveplots:
    plt.savefig(os.path.join(output_dir_base, "images.png"), facecolor="white", bbox_inches="tight", dpi=150)


# 保存数据到本地
if save_data:
    # 1. 保存观测表
    cadc_tab.write(os.path.join(output_dir_base, f"observations_{source}.csv"), format='ascii.csv', overwrite=True)
    cadc_tab.write(os.path.join(output_dir_base, f"observations_{source}.fits"), format='fits', overwrite=True)
    print(f"✅ 观测数据表已保存至: {output_dir_base}/")

    # 2. 保存检测结果表（如果存在）
    if dtab is not None and len(dtab) > 0:
        dtab.write(os.path.join(output_dir_base, f"detections_{source}.csv"), format='ascii.csv', overwrite=True)
        dtab.write(os.path.join(output_dir_base, f"detections_{source}.fits"), format='fits', overwrite=True)
        print(f"✅ 检测结果表已保存至: {output_dir_base}/")

    # 3. 保存 PNG 缩略图
    png_dir = os.path.join(output_dir_base, "thumbnails")
    os.makedirs(png_dir, exist_ok=True)
    for i, img in enumerate(goodimages):
        mjd = goodtab['MJD'][i]
        filt = goodtab['Filter'][i]
        path = os.path.join(png_dir, f"img_{i:03d}_{filt}_{mjd:.5f}.png")
        img.save(path)
    print(f"✅ {ngood} 张缩略图已保存至: {png_dir}/")

    # 4. 下载原始 warp FITS 文件
    fits_dir = os.path.join(output_dir_base, "fits_data")
    print(f"\n开始下载 {len(goodtab)} 个原始 warp FITS 文件...")
    downloaded_fits = download_warp_fits(goodtab, fits_dir)
    print(f"✅ 共下载 {len(downloaded_fits)} 个 FITS 文件到: {fits_dir}/")


print("\n" + "=" * 60)
print("If you use CADC facilities or SSOS in your research, please include:")
print("This research used the facilities of the Canadian Astronomy Data Centre")
print("operated by the National Research Council of Canada with the support")
print("of the Canadian Space Agency.")
print("Please cite: Gwyn, Hill & Kavelaars (2012)")
print("=" * 60)