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
from tkinter import simpledialog


# 配置选项
saveplots = True  # 是否保存生成的图表
save_data = True  # 是否将观测表和图像保存到本地
use_jpl = True  # 是否尝试从 JPL Horizons 获取高精度星历（设为 False 跳过）
source = "2016 OR17"  # 目标编号

# 网络请求/下载配置
DOWNLOAD_MAX_ATTEMPTS = 5  # 包含第一次请求在内的总尝试次数
CONNECT_TIMEOUT = 15  # 建立连接超时（秒）
READ_TIMEOUT = 120  # 数据读取超时（秒）
DOWNLOAD_RETRY_BACKOFF = 2.0  # 指数退避基数：2、4、8、16 ... 秒
RETRYABLE_HTTP_STATUS = {408, 429, 500, 502, 503, 504}

astropy.conf.max_width = 150


def _retry_wait(attempt, backoff=DOWNLOAD_RETRY_BACKOFF):
    """返回第 attempt 次失败后的指数退避等待秒数。"""
    return backoff * (2 ** (attempt - 1))


def request_get_with_retry(url, *, max_attempts=DOWNLOAD_MAX_ATTEMPTS,
                           timeout=(CONNECT_TIMEOUT, READ_TIMEOUT),
                           backoff=DOWNLOAD_RETRY_BACKOFF, **kwargs):
    """执行带重试的 GET 请求。

    连接错误、超时以及常见临时 HTTP 错误会自动重试；普通 4xx（例如 404）
    会立即抛出，避免对永久性错误进行无意义重试。
    """
    last_error = None

    for attempt in range(1, max_attempts + 1):
        try:
            response = requests.get(url, timeout=timeout, **kwargs)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as exc:
            last_error = exc
            status = exc.response.status_code if exc.response is not None else None
            retryable = status is None or status in RETRYABLE_HTTP_STATUS

            if not retryable or attempt >= max_attempts:
                raise

            wait = _retry_wait(attempt, backoff)
            print(f"\n请求失败（第 {attempt}/{max_attempts} 次）: {exc}")
            print(f"{wait:.1f} 秒后重试: {url}")
            time.sleep(wait)

    raise last_error


def init_casjobs():
    """初始化可选的 MAST CasJobs 连接。

    未提供凭据时仍保留原来的图形输入方式；用户取消输入或连接失败时不再终止
    主程序，而是跳过后续目录匹配。
    """
    if not os.environ.get('CASJOBS_USERID') or not os.environ.get('CASJOBS_PW'):
        root = None
        try:
            root = tk.Tk()
            root.withdraw()

            userid = simpledialog.askstring("CasJobs Login", "Enter Casjobs UserID:")
            password = None
            if userid and userid.strip():
                password = simpledialog.askstring("CasJobs Login", "Enter Casjobs password:", show='*')

            if not userid or not userid.strip() or not password or not password.strip():
                print("Warning: CasJobs credentials were not provided. Skipping CasJobs catalog queries.")
                return None

            os.environ['CASJOBS_USERID'] = userid.strip()
            os.environ['CASJOBS_PW'] = password.strip()
        except Exception as exc:
            print(f"Warning: Could not request CasJobs credentials: {exc}")
            print("Continuing without CasJobs catalog queries.")
            return None
        finally:
            if root is not None:
                try:
                    root.destroy()
                except tk.TclError:
                    pass

    try:
        jobs = mastcasjobs.MastCasJobs(context="MyDB")
        print("Confirmed MAST CasJobs account access")
        return jobs
    except Exception as exc:
        print(f"Warning: Failed to connect to CasJobs: {exc}")
        print("Continuing without CasJobs catalog queries.")
        return None


def cadc_ssos_query(object_name, search="bynameCADC", epoch1=54985, epoch2=57079,
                    xyres="no", telinst="Pan-STARRS1", lang="en", format="tsv",
                    url="https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl"):
    """使用 CADC 移动天体查询功能按目标名称查找 PS1 观测数据。"""
    params = dict(lang=lang, object=object_name, search=search,
                  epoch1=epoch1, epoch2=epoch2, xyres=xyres,
                  telinst=telinst, format=format)

    with request_get_with_retry(url, params=params) as response:
        text = response.text

    msg = text.split('\n', 1)[0]
    if "not found" in msg.lower():
        if msg.endswith("<br/>"):
            msg = msg[:-5]
        raise ValueError(f"Error from CADC search: {msg}")

    return Table.read(text, format="ascii.csv", delimiter="\t",
                      converters=dict(Image_target=str))


def getimages(tab, format="jpg", verbose=False, **kw):
    """下载表格指定的 PS1 cutout 图像。

    单张图像的临时网络错误会自动重试。达到最大尝试次数后，将该图像记为
    None 并继续处理后续图像，避免一次失败导致整个脚本退出。
    """
    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be jpg, png or fits")

    urls = tab2cutout(tab, format=format, verbose=verbose, **kw)
    images = []

    if verbose:
        t0 = time.time()
        print(f"{time.time() - t0:.1f} s: completed 0 of {len(urls)} images", end="")

    for i, url in enumerate(urls, start=1):
        try:
            with request_get_with_retry(url) as response:
                content = response.content

            if format == "fits":
                images.append(fits.open(BytesIO(content), memmap=False))
            else:
                img = Image.open(BytesIO(content))
                img.load()
                images.append(img.copy())
                img.close()

        except (requests.exceptions.RequestException, UnidentifiedImageError, OSError) as exc:
            print(f"\n图像下载失败，已跳过: {url}")
            print(f"原因: {exc}")
            images.append(None)

        if verbose:
            print(f"\r{time.time() - t0:.1f} s: completed {i} of {len(urls)} images", end="")

    if verbose:
        print()

    return images


def tab2path(tab, imagecol='Image', verbose=False):
    """将包含常规列的表格转换为 warp FITS 路径。"""
    for c in ('projcell', 'subcell', imagecol):
        assert c in tab.colnames, f"'{c}' column not found in table"

    rv = []
    for row in tab:
        projcell = f"{row['projcell']:04d}"
        subcell = f"{row['subcell']:03d}"
        image = row[imagecol]
        rv.append(f"/rings.v3.skycell/{projcell}/{subcell}/{image}.fits")
    return rv


def tab2cutout(tab, size=64, output_size=256, format="jpg", autoscale=99.5, asinh=True,
               racol="ra", deccol="dec", imagecol="Image", verbose=False,
               urlbase="https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"):
    """将包含常规列的表格转换为 warp 图像裁剪 URL。"""
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
    d1 = np.minimum(np.minimum(x1 + nx1 // 2, nx1 // 2 - 1 - x1),
                    np.minimum(y1 + ny1 // 2, ny1 // 2 - 1 - y1))
    d2 = np.minimum(np.minimum(x2 + nx2 // 2, nx2 // 2 - 1 - x2),
                    np.minimum(y2 + ny2 // 2, ny2 // 2 - 1 - y2))
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


def download_warp_fits(tab, download_dir,
                       base_url="https://ps1images.stsci.edu/data/ps1/data/",
                       max_attempts=DOWNLOAD_MAX_ATTEMPTS,
                       timeout=(CONNECT_TIMEOUT, READ_TIMEOUT),
                       backoff=DOWNLOAD_RETRY_BACKOFF):
    """下载原始 warp FITS 文件，失败时自动重试并继续后续文件。

    下载先写入 ``.part`` 临时文件，完整写入并通过 Content-Length（若服务器
    提供）校验后再原子改名，从而避免残缺文件在下次运行时被误认为已完成。
    """
    os.makedirs(download_dir, exist_ok=True)
    paths = tab2path(tab, imagecol='Image')
    downloaded_files = []
    failed_files = []

    for i, path in enumerate(paths):
        filename = os.path.basename(path)
        url = f"{base_url.rstrip('/')}/{path.lstrip('/')}"
        filepath = os.path.join(download_dir, filename)
        part_path = filepath + ".part"

        if os.path.exists(filepath) and os.path.getsize(filepath) > 0:
            print(f"  已存在，跳过: {filepath}")
            downloaded_files.append(filepath)
            continue

        if os.path.exists(filepath) and os.path.getsize(filepath) == 0:
            os.remove(filepath)

        print(f"[{i + 1}/{len(tab)}] 下载: {url}")
        success = False

        for attempt in range(1, max_attempts + 1):
            try:
                if os.path.exists(part_path):
                    os.remove(part_path)

                with requests.get(url, stream=True, timeout=timeout) as response:
                    response.raise_for_status()
                    expected_size = response.headers.get('Content-Length')
                    expected_size = int(expected_size) if expected_size and expected_size.isdigit() else None

                    bytes_written = 0
                    with open(part_path, 'wb') as fh:
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            if not chunk:
                                continue
                            fh.write(chunk)
                            bytes_written += len(chunk)

                if bytes_written == 0:
                    raise OSError("服务器返回了空文件")
                if expected_size is not None and bytes_written != expected_size:
                    raise OSError(
                        f"文件大小校验失败：预期 {expected_size} 字节，实际 {bytes_written} 字节"
                    )

                os.replace(part_path, filepath)
                print(f" 成功保存: {filepath}")
                downloaded_files.append(filepath)
                success = True
                break

            except requests.exceptions.RequestException as exc:
                if os.path.exists(part_path):
                    os.remove(part_path)

                status = exc.response.status_code if exc.response is not None else None
                retryable = status is None or status in RETRYABLE_HTTP_STATUS

                if not retryable or attempt >= max_attempts:
                    print(f" 下载失败，已放弃 {url}: {exc}")
                    break

                wait = _retry_wait(attempt, backoff)
                print(f" 下载失败（第 {attempt}/{max_attempts} 次）: {exc}")
                print(f" {wait:.1f} 秒后重试...")
                time.sleep(wait)

            except OSError as exc:
                if os.path.exists(part_path):
                    os.remove(part_path)

                if attempt >= max_attempts:
                    print(f" 保存/校验失败，已放弃 {filepath}: {exc}")
                    break

                wait = _retry_wait(attempt, backoff)
                print(f" 保存/校验失败（第 {attempt}/{max_attempts} 次）: {exc}")
                print(f" {wait:.1f} 秒后重试...")
                time.sleep(wait)

        if not success:
            failed_files.append(url)

    if failed_files:
        print(f"\n⚠️ 有 {len(failed_files)} 个 FITS 文件最终下载失败；其他文件已继续处理。")
        for failed_url in failed_files:
            print(f"  - {failed_url}")

    return downloaded_files


def main():
    jobs = init_casjobs()

    output_dir_base = f"PS1_data_{source}"
    os.makedirs(output_dir_base, exist_ok=True)

    cadc_tab = cadc_ssos_query(source)
    print(f"Returned table has {len(cadc_tab)} rows")
    print(cadc_tab[:8])

    ra_col = cadc_tab['Object_RA'].copy()
    dec_col = cadc_tab['Object_Dec'].copy()
    fullname = f"Object {source}"

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
        except Exception as exc:
            print(f"Warning: Could not retrieve data from JPL Horizons: {exc}")
            print("Falling back to CADC-provided positions.")

    for c in ('ra', 'dec'):
        values = [ra_col, dec_col][['ra', 'dec'].index(c)]
        if c in cadc_tab.colnames:
            cadc_tab[c] = values
        else:
            cadc_tab.add_column(values, name=c, index=cadc_tab.colnames.index('Object_RA'))

    cadc_tab['dra'] = (cadc_tab['ra'] - cadc_tab['Object_RA']) * 3600 * np.cos(np.radians(cadc_tab['dec']))
    cadc_tab['ddec'] = (cadc_tab['dec'] - cadc_tab['Object_Dec']) * 3600
    for c in ('ra', 'dec', 'Object_RA', 'Object_Dec'):
        cadc_tab[c].format = ".6f"
    for c in ('dra', 'ddec'):
        cadc_tab[c].format = ".2f"
        cadc_tab[c].unit = "arcsec"

    print(f"dra range  {cadc_tab['dra'].min():.2f} {cadc_tab['dra'].max():.2f} "
          f"MAD {np.median(np.abs(cadc_tab['dra'])):.2f} arcsec")
    print(f"ddec range {cadc_tab['ddec'].min():.2f} {cadc_tab['ddec'].max():.2f} "
          f"MAD {np.median(np.abs(cadc_tab['ddec'])):.2f} arcsec")

    rv = findskycell(cadc_tab['ra'], cadc_tab['dec'])
    for c in rv.colnames:
        if c not in cadc_tab.colnames:
            cadc_tab[c] = rv[c]
    cadc_tab["crval1"].format = ".6f"
    for c in ("x", "y"):
        cadc_tab[c].format = ".2f"

    umjd, _ = np.unique(cadc_tab['MJD'], return_counts=True)
    nuniq = len(umjd)
    print(f"cadc_tab has {len(cadc_tab) - nuniq} duplicate entries, leaving {nuniq} unique entries")

    keep = np.array([
        image.startswith(f"rings.v3.skycell.{projcell:04d}.{subcell:03d}.")
        for image, projcell, subcell in
        zip(cadc_tab['Image'].data, cadc_tab['projcell'].data, cadc_tab['subcell'].data)
    ])

    if keep.sum() != nuniq:
        jtest = join(cadc_tab['MJD', 'Image', 'projcell', 'subcell'],
                     cadc_tab['MJD', 'Image_target'][keep], join_type='left')
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
                keep = np.array([
                    image.startswith(f"rings.v3.skycell.{projcell:04d}.{subcell:03d}.")
                    for image, projcell, subcell in
                    zip(cadc_tab['Image'].data, cadc_tab['projcell'].data, cadc_tab['subcell'].data)
                ])
                jtest = join(cadc_tab['MJD', 'Image', 'projcell', 'subcell'],
                             cadc_tab['MJD', 'Image_target'][keep], join_type='left')
                w = np.where(jtest['Image_target'].mask)[0]
        assert len(w) == nuniq - keep.sum()
        assert (cadc_tab['Image'][w] == jtest['Image'][w]).all()
        new_projcell = [int(x.split('.')[0]) for x in cadc_tab['Image_target'][w].data]
        new_subcell = [int(x.split('.')[1]) for x in cadc_tab['Image_target'][w].data]
        print(f"Patching projcell and subcell for {len(w)} images")
        cadc_tab['projcell'][w] = new_projcell
        cadc_tab['subcell'][w] = new_subcell
        keep = np.array([
            image.startswith(f"rings.v3.skycell.{projcell:04d}.{subcell:03d}.")
            for image, projcell, subcell in
            zip(cadc_tab['Image'].data, cadc_tab['projcell'].data, cadc_tab['subcell'].data)
        ])
        assert keep.sum() == nuniq
        print("Modified lines:")
        print(cadc_tab['Image', 'MJD', 'Filter', 'Image_target', 'projcell', 'subcell'][w])

    cadc_tab = cadc_tab[keep]

    emin = Time(54985.0, format='mjd')
    emax = Time(57079.0, format='mjd')
    print(f"Searching for {source} from {emin.isot} to {emax.isot} ({(emax.mjd - emin.mjd):.0f} days)")

    plt.rcParams.update({"font.size": 14})
    plt.figure(1, (12, 6))

    ratab = cadc_tab['ra']
    dec_plot = cadc_tab['dec']
    jpl_mjd = cadc_tab['MJD']

    if use_jpl:
        try:
            obj = Horizons(id=source, location='F51',
                           epochs=dict(start=emin.isot, stop=emax.isot, step='1d'))
            full_jpltab = obj.ephemerides(extra_precision=True)
            jpl_mjd = Time(full_jpltab['datetime_jd'], format='jd').mjd
            ratab = full_jpltab['RA'].data.copy()
            dec_plot = full_jpltab['DEC'].data.copy()
            w = np.where(ratab[1:] > ratab[:-1] + 180)[0]
            for k in w:
                ratab[k + 1:] -= 360
            w = np.where(ratab[1:] < ratab[:-1] - 180)[0]
            for k in w:
                ratab[k + 1:] += 360
            while ratab.min() < 0:
                ratab += 360
            plt.plot(ratab, dec_plot, '-', alpha=0.6, color='gray', label="Orbital path")
            mstep = max(1, int(30 / (jpl_mjd[1] - jpl_mjd[0]) + 0.5))
            plt.plot(ratab[::mstep], dec_plot[::mstep], 'o', color="tab:blue", label="Monthly positions")
        except Exception as exc:
            print(f"Could not get full ephemeris: {exc}")

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

    imsize = 65
    output_size = 130
    images = getimages(cadc_tab, size=imsize, output_size=output_size, verbose=True)
    print(f"Got {len(images)} image slots")

    igood = [
        i for i, img in enumerate(images)
        if img is not None and np.array(img).min() != 255
    ]
    goodimages = [images[i] for i in igood]
    ngood = len(goodimages)
    goodtab = cadc_tab[igood]
    goodtab.add_column(np.arange(len(goodtab)), name="SrcID", index=0)
    print(f"{ngood} of {len(images)} images downloaded successfully and are not blank")

    dtab = None
    if jobs is not None and len(goodtab) > 0:
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
                plt.savefig(os.path.join(output_dir_base, "position-diffs.png"),
                            facecolor="white", bbox_inches="tight")

        except Exception as exc:
            print(f"Warning: Could not query PS1 catalog: {exc}")
            print("Continuing without source position markings.")
            dtab = None

    if ngood > 0:
        ncols = 7
        nrows = (ngood + ncols - 1) // ncols
        hsize = ncols * 2.5
        vsize = nrows * 2.5
        extent = np.array([1, -1, -1, 1]) * 0.25 * imsize / 2

        plt.rcParams.update({"font.size": 10})
        plt.figure(2, (hsize, vsize))

        goodmatch = None
        if dtab is not None and len(dtab) > 0:
            goodmatch = join(goodtab, dtab['SrcID', 'ra', 'dec', 'psfFlux', 'psfQfPerfect'],
                             keys='SrcID', join_type='left')
            mask = ~np.isnan(goodmatch['psfFlux'])
            mag = np.full(len(goodmatch), np.nan)
            mag[mask] = -2.5 * np.log10(goodmatch['psfFlux'][mask]) + 8.90
            goodmatch['mag'] = mag

        for i in range(ngood):
            plt.subplot(nrows, ncols, i + 1)
            plt.imshow(goodimages[i], origin="upper", cmap="gray", vmin=0, vmax=255, extent=extent)

            if (goodmatch is not None and i < len(goodmatch)
                    and not np.isnan(goodmatch['ra_2'][i])
                    and not np.isnan(goodmatch['dec_2'][i])):
                dx_val = (goodmatch['ra_2'][i] - goodmatch['ra_1'][i]) * 3600 * np.cos(
                    np.radians(goodmatch['dec_1'][i]))
                dy_val = (goodmatch['dec_2'][i] - goodmatch['dec_1'][i]) * 3600
                color = "lightblue" if goodmatch['psfQfPerfect'][i] > 0.95 else "red"
                plt.plot(dx_val, dy_val, 'o', color=color, markerfacecolor="none",
                         markersize=15, markeredgewidth=2)
                if not np.isnan(goodmatch['mag'][i]):
                    plt.text(dx_val, dy_val - 5, f"{goodmatch['mag'][i]:.2f}", color="yellow",
                             fontsize=8, ha='center', va='top', weight='bold')

            plt.title(f"{goodtab['Filter'][i]} {goodtab['MJD'][i]:.5f}")

        plt.tight_layout()
        plt.suptitle(
            f"{ngood} {imsize}x{imsize} pixel ({imsize / 4:.0f}x{imsize / 4:.0f} arcsec) "
            f"PS1 images for {fullname}", y=1.01, fontsize=16, fontweight="bold"
        )
        if saveplots:
            plt.savefig(os.path.join(output_dir_base, "images.png"),
                        facecolor="white", bbox_inches="tight", dpi=150)
    else:
        print("Warning: No usable thumbnail images were downloaded; skipping image grid.")

    if save_data:
        cadc_tab.write(os.path.join(output_dir_base, f"observations_{source}.csv"),
                       format='ascii.csv', overwrite=True)
        cadc_tab.write(os.path.join(output_dir_base, f"observations_{source}.fits"),
                       format='fits', overwrite=True)
        print(f"✅ 观测数据表已保存至: {output_dir_base}/")

        if dtab is not None and len(dtab) > 0:
            dtab.write(os.path.join(output_dir_base, f"detections_{source}.csv"),
                       format='ascii.csv', overwrite=True)
            dtab.write(os.path.join(output_dir_base, f"detections_{source}.fits"),
                       format='fits', overwrite=True)
            print(f"✅ 检测结果表已保存至: {output_dir_base}/")

        png_dir = os.path.join(output_dir_base, "thumbnails")
        os.makedirs(png_dir, exist_ok=True)
        for i, img in enumerate(goodimages):
            mjd = goodtab['MJD'][i]
            filt = goodtab['Filter'][i]
            path = os.path.join(png_dir, f"img_{i:03d}_{filt}_{mjd:.5f}.png")
            img.save(path)
        print(f"✅ {ngood} 张缩略图已保存至: {png_dir}/")

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


if __name__ == "__main__":
    main()
