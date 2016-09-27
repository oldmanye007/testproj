"""
test for documentation
"""

import sys, os, re, math, struct, argparse, csv
import envi_header_handler as EHH #Custom program written by Steve Cochaver
import numpy as np

try:
    import gdal
except:
    try:
        from osgeo import gdal
    except:
        print "Import of gdal failed."
        sys.exit(1)

from gdalconst import *


def gdaltype2nptype(gdalDataType):
  """
     test for gdaltype2nptype
  """
  typelist=['B','H','h','I','i','f']
  fmt=typelist[gdalDataType-1]
  return fmt

def remove_hdr_item(hdr, itemkey):
  """
    test for remove_hdr_item
  """
  if itemkey in hdr._key_order: hdr._key_order.remove(itemkey)
  if itemkey in hdr._nested_keys: hdr._nested_keys.remove(itemkey)
  hdr.hdr_dict.pop(itemkey, None)
  
def getraster(band):
  """
    test for getraster
  """
  typelist=['B','H','h','I','i','f'] 
  scanline = band.ReadRaster( 0, 0, band.XSize, band.YSize, band.XSize, band.YSize, band.DataType )
  fmt=typelist[band.DataType-1]

  value = struct.unpack(fmt * band.XSize*band.YSize, scanline)
  # http://geoinformaticstutorial.blogspot.com/2012/09/reading-raster-data-with-python-and-gdal.html
  return np.array(value)  
  
def find_col(target, header):
  index=0
  for item in header:
    if target == item:
      return index
    index+=1
  if index == len(header):
    return -1


def read_csv_coef(incsv):

  #csv_matrix=np.loadtxt(incsv,delimiter=",",skiprows=1,usecols=(1,2))
  #print csv_matrix.shape
  #sys.exit(1)
  with open(incsv, 'rb') as csvfile:
    r = csv.reader(csvfile)
    #print csvfile.fieldnames
    header = r.next()
    print 'csv header',header,
    num_col=len(header)
    r.next()
    #csvfile.seek(0)

    col_tuple=tuple(range(1,num_col))

    #sys.exit(1)
    csv_matrix=np.loadtxt(incsv,skiprows=1,delimiter=",",usecols=col_tuple)
    #csv_matrix=np.loadtxt(incsv,skiprows=1, usecols=col_tuple)
    #print csv_matrix.shape
    
    band_list = [int(row[0][5:]) for row in r]
    #print type(band_list)
    #print band_list[0:4]
    
    return {"spec_name":header[1:],"data":csv_matrix,"bandlist":band_list}

    
def cal_pls_y(inimg,outdir, csvdata):

  v_nodata=-50

  inBN=os.path.basename(os.path.splitext(inimg)[0])
  outimg=outdir+'/'+inBN+'_QUERPRN_4models'
  print outimg
  
  if os.path.exists(inimg+".hdr"):
    hdrinfo=EHH.ENVI_Header(inimg+".hdr")
    print inimg+".hdr"
  else:
    if os.path.exists(os.path.splitext(inimg)[0]+".hdr"):
      hdrinfo=EHH.ENVI_Header(os.path.splitext(inimg)[0]+".hdr")
      print os.path.splitext(inimg)[0]+".hdr"
    else:
      print "Cannot find .hdr file"
      sys.exit(1)
  
  nband=int(hdrinfo.get_value('bands'))
  m_interleave=hdrinfo.get_value('interleave')
  
  print len(csvdata["spec_name"]),' species'
  
  in_ds=gdal.Open(inimg, GA_ReadOnly)

  if in_ds is None:
    print 'Could not open ' + inimg
    sys.exit(1)
  try: 
     
    if (os.path.exists(outimg)):
      print "File exists!"
      in_ds=None   
      sys.exit(1)

    nband=in_ds.RasterCount
    print nband, in_ds.RasterYSize, in_ds.RasterXSize  

    remove_hdr_item(hdrinfo,'wavelength units' )
    remove_hdr_item(hdrinfo,'wavelength' )
    remove_hdr_item(hdrinfo,'bbl' )
    remove_hdr_item(hdrinfo,'z plot titles' )
    remove_hdr_item(hdrinfo,'z plot range' )
    remove_hdr_item(hdrinfo,'fwhm' )
    remove_hdr_item(hdrinfo,'default bands' )

    csv_coef=csvdata["data"]
    intercept_list=np.reshape(csv_coef[0,],(len(csvdata["spec_name"])))
    coef_matrix=csv_coef[1:,]
    total_band=coef_matrix.shape[0]
    coef_matrix=np.reshape(coef_matrix,(total_band,len(csvdata["spec_name"])))
    
    
    print coef_matrix.shape
    print intercept_list
    #sys.exit(1)
    with open(outimg, 'ab') as f:

      if m_interleave=='bsq':
        print m_interleave    
        i_spec=0
        for species in csvdata["spec_name"]:
          print species
          
          #sumband=np.ones((in_ds.RasterYSize * in_ds.RasterXSize),np.float)
          sumband=np.zeros((in_ds.RasterYSize * in_ds.RasterXSize),np.float)
          #print csvdata
          #print intercept_list[i_spec]
          #print sumband[0:10]
          
          #print sumband.shape, sumband[0:10]
          #sys.exit(1)
          vec_norm=np.zeros((in_ds.RasterYSize * in_ds.RasterXSize),np.float)
          nodata_mask=np.zeros((in_ds.RasterYSize * in_ds.RasterXSize),np.float)
          #sum_coef=0
          
          for index in range(total_band):
            iband=csvdata["bandlist"][index]
            #print index, iband,i_spec
            cur_coef= coef_matrix[index,i_spec]
            band=in_ds.GetRasterBand(iband)
            imgband=getraster(band)
            tmp=imgband*cur_coef
            #sum_coef+=cur_coef
            #if index < 10:             
            #print cur_coef, index, iband
            vec_norm=np.add(vec_norm,np.square(imgband))
            sumband=np.add(sumband,tmp)        
            
            if (index== 0):
              nodata_mask[imgband ==v_nodata]=1
            #sys.exit(1)
            
          #v_nodata=v_nodata*sum_coef+ intercept_list[i_spec]
          #v_nodata= sum_coef*v_nodata/abs(v_nodata)/math.sqrt(total_band)+intercept_list[i_spec]
          vec_norm=np.sqrt(vec_norm)

          sumband=np.divide(sumband,vec_norm)
          #print sumband[0],vec_norm[0],intercept_list[i_spec]
          sumband=sumband+intercept_list[i_spec]
          sumband[nodata_mask==1]=-9999
          #sumband[np.where(sumband<0) and np.where(sumband>-9999)]=0
          f.write(sumband.astype(np.float32))
          #print v_nodata
          i_spec+=1
          #break
      hdrinfo.change_value("interleave", "bsq")  
      hdrinfo.change_value("bands", str(len(csvdata["spec_name"])))
      hdrinfo.change_value("data type", '4')
      hdrinfo.change_value("band names",csvdata["spec_name"])    
      hdrinfo.write_header(outdir, outimg+".hdr")
       
  except RuntimeError, e:
    print e
    in_ds=None 
    sys.exit(1)  
  
  in_ds=None       
    

def main(argv):
# python  h:/plot/apply_pls.py -i "h:/plot/tmp_serak/f090714t01p00r05rdn_b_ort_img_tafkaa_orig_refl_img_bsq_trc_xtr" -c h:/plot/rawest_test.csv -o h:/plot/out/

# python  h:/plot/apply_pls.py -i "h:/plot/tmp_serak/f090714t01p00r06rdn_b_ort_img_tafkaa_orig_refl_img_bsq_trc_xtr" -c h:/plot/rawest_test.csv -o h:/plot/out/

# python  h:/plot/apply_pls.py -i "h:/plot/tmp_serak/f090714t01p00r05rdn_b_ort_img_tafkaa_orig_refl_img_bsq_trc_xtr" -c h:/plot/pidQUERPRN_PLS0_4models.csv -o h:/plot/out/

# python  h:/plot/apply_pls.py -i "h:/plot/tmp_serak/f090714t01p00r06rdn_b_ort_img_tafkaa_orig_refl_img_bsq_trc_xtr" -c h:/plot/pidQUERPRN_PLS0_4models.csv -o h:/plot/out/


  parser = argparse.ArgumentParser(description='This code is for PLS coefficients applying to AVIRIS image')
  parser.add_argument('-i','--inimg', type=str, help='Input multi-band image file name',required=True)
  parser.add_argument('-c','--coefcsv', type=str, help='Input csv file of coefficients.',required=True)
  parser.add_argument('-o','--outdir', type=str, help='Output dir', required=True)  #OBS_ORT image

  args = parser.parse_args()
  
  coefcsv=args.coefcsv
  inimg=args.inimg  # default bsq 224 bands
  outdir=args.outdir
  
  csvdata=read_csv_coef(coefcsv)  # include intercepts
  print csvdata["spec_name"]
  print csvdata["data"].shape
  print csvdata["bandlist"][0:10]  #band number starts from 1

  cal_pls_y(inimg,outdir, csvdata)


if __name__ == "__main__":
  """
    Main function
  """
  main(sys.argv[1:]) 