from mimetypes import init
import os
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.spatial import cKDTree
import shapefile as shp
import datetime
from flopy.utils.gridintersect import GridIntersect
from flopy.utils import Raster
from shapely.geometry import Polygon, Point, LineString
import rasterio
from rasterio.features import rasterize
import geopandas as gpd
from flopy.mf6.modflow.mfgwfnpf import ModflowGwfnpf
from flopy.mf6.modflow.mfgwfsto import ModflowGwfsto
pd.options.mode.chained_assignment = None  # default='warn'



def set_hydraulic_properties(sim, gdf, s, bbox, delr, delc, hk_attributes, vk_attributes, sy_attributes, ss_attributes):
    # Define the spatial resolution of your MODFLOW grid
    xres = delr
    yres = delc

    # Define the bounds of your MODFLOW grid
    xmin = int((bbox[0])//100*100)
    ymin = int((bbox[1])//100*100+100)
    xmax = int((bbox[2])//100*100)
    ymax = int((bbox[3])//100*100+100)

    # Create an empty raster dataset with the same spatial resolution and bounds as your MODFLOW grid
    transform = rasterio.transform.from_origin(xmin, ymax, xres, yres)
    out_shape = (int((ymax-ymin)/yres), int((xmax-xmin)/xres))

    # Loop over each layer
    for i in range(len(hk_attributes)):
        # Check if there are any non-zero values in the shapefile for this attribute
        if any(gdf[hk_attributes[i]] != 0):
            # Rasterize the horizontal hydraulic conductivity attribute of the shapefile for this layer
            hk_raster = rasterize(
                ((geom,value) for geom, value in zip(gdf.geometry, gdf[hk_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,  # use 1 as the fill value
                default_value=1,
                dtype='float64'
            )
        else:
            # If all values are 0, create a raster filled with the default value
            hk_raster = np.full(out_shape, 1, dtype='float64')

        # Set the horizontal hydraulic conductivity in this layer of the MODFLOW model
        sim.lpf.hk[i] = hk_raster

        # Repeat similar steps for 'vka', 'sy', and 'ss'
        if any(gdf[vk_attributes[i]] != 0):
            vk_raster = rasterize(
                ((geom,value) for geom, value in zip(gdf.geometry, gdf[vk_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,
                default_value=0.1,
                dtype='float64'
            )
        else:
            vk_raster = np.full(out_shape, 0.1, dtype='float64')

        sim.lpf.vka[i] = vk_raster

        if any(s[sy_attributes[i]] != 0):
            sy_raster = rasterize(
                ((geom,value) for geom, value in zip(s.geometry, s[sy_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,
                default_value=0.05,
                dtype='float64'
            )
        else:
            sy_raster = np.full(out_shape, 0.05, dtype='float64')

        sim.lpf.sy[i] = sy_raster

        if any(s[ss_attributes[i]] != 0):
            ss_raster = rasterize(
                ((geom,value) for geom, value in zip(s.geometry, s[ss_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,
                default_value=0.001,
                dtype='float64'
            )
        else:
            ss_raster = np.full(out_shape, 0.001, dtype='float64')

        # Check if any values in the raster exceed the threshold
        ss_raster = np.where(ss_raster > 0.009, 0.009, ss_raster)

        # Set the specific storage in this layer of the MODFLOW model
        sim.lpf.ss[i] = ss_raster

from flopy.mf6.modflow.mfgwfnpf import ModflowGwfnpf

from flopy.mf6.modflow.mfgwfnpf import ModflowGwfnpf
from flopy.mf6.modflow.mfgwfsto import ModflowGwfsto

def set_hydraulic_properties_mf6(sim, gwf, gdf, s, bbox, delr, delc, hk_attributes, vk_attributes, sy_attributes, ss_attributes):
    # Define the spatial resolution of your MODFLOW grid
    xres = delr
    yres = delc

    # Define the bounds of your MODFLOW grid
    xmin = int((bbox[0])//100*100)
    ymin = int((bbox[1])//100*100+100)
    xmax = int((bbox[2])//100*100)
    ymax = int((bbox[3])//100*100+100)

    # Create an empty raster dataset with the same spatial resolution and bounds as your MODFLOW grid
    transform = rasterio.transform.from_origin(xmin, ymax, xres, yres)
    out_shape = (int((ymax-ymin)/yres), int((xmax-xmin)/xres))

    # Get the number of layers, rows, and columns
    nlay, nrow, ncol = gwf.dis.nlay.array, gwf.dis.nrow.array, gwf.dis.ncol.array

    # Initialize k, k33, sy, and ss as arrays
    k_array = np.ones((nlay, nrow, ncol))
    k33_array = np.ones((nlay, nrow, ncol))
    sy_array = np.ones((nlay, nrow, ncol))
    ss_array = np.ones((nlay, nrow, ncol))

    # Loop over each layer
    for i in range(len(hk_attributes)):
        # Check if there are any non-zero values in the shapefile for this attribute
        if any(gdf[hk_attributes[i]] != 0):
            # Rasterize the horizontal hydraulic conductivity attribute of the shapefile for this layer
            hk_raster = rasterize(
                ((geom,value) for geom, value in zip(gdf.geometry, gdf[hk_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,  # use 1 as the fill value
                default_value=1,
                dtype='float64'
            )
        else:
            # If all values are 0, create a raster filled with the default value
            hk_raster = np.full(out_shape, 1, dtype='float64')

        # Set the horizontal hydraulic conductivity in this layer of the MODFLOW model
        k_array[i] = hk_raster

        # Repeat similar steps for 'vka', 'sy', and 'ss'
        if any(gdf[vk_attributes[i]] != 0):
            vk_raster = rasterize(
                ((geom,value) for geom, value in zip(gdf.geometry, gdf[vk_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,
                default_value=0.1,
                dtype='float64'
            )
        else:
            vk_raster = np.full(out_shape, 0.1, dtype='float64')

        k33_array[i] = vk_raster

        if any(s[sy_attributes[i]] != 0):
            sy_raster = rasterize(
                ((geom,value) for geom, value in zip(s.geometry, s[sy_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,
                default_value=0.05,
                dtype='float64'
            )
        else:
            sy_raster = np.full(out_shape, 0.05, dtype='float64')

        sy_array[i] = sy_raster

        if any(s[ss_attributes[i]] != 0):
            ss_raster = rasterize(
                ((geom,value) for geom, value in zip(s.geometry, s[ss_attributes[i]]) if value != 0),
                out_shape=out_shape,
                transform=transform,
                fill=1,
                default_value=0.001,
                dtype='float64'
            )
        else:
            ss_raster = np.full(out_shape, 0.001, dtype='float64')

        # Check if any values in the raster exceed the threshold
        ss_raster = np.where(ss_raster > 0.009, 0.009, ss_raster)

        # Set the specific storage in this layer of the MODFLOW model
        ss_array[i] = ss_raster

    # Create a new npf package with the k and k33 values
    npf = ModflowGwfnpf(gwf, k=k_array, k33=k33_array)

    # Create a new sto package with the sy and ss values
    sto = ModflowGwfsto(gwf, sy=sy_array, ss=ss_array, steady_state={0: True}, save_flows=False)

def interpolate_to_grid(file_list,xul,yul,nrow,ncol,delr,delc,skiprows=None):
    """interpolate asc data to a model grid using natural neighbor

    Args:
        file_list ('list'): list of input files 
        nrows ('int'): model number of rows
        ncols ('int'): model number of rows
        delr ('float'): model delta row
        delc ('float'): model delta column
        skiprows('int'): number of rows to skip in file

    Return:
        3D array with data interpolated
    """

    # a model grid of (x,y) coordinates is generated
    model_xy = np.empty((nrow*ncol, 2), dtype=float)

    inode = 0
    for irow in range(nrow):
        for icol in range(ncol):
            model_xy[inode,0]=xul+delc*(icol+1/2)
            model_xy[inode,1]=yul-delr*(irow+1/2)
            inode=inode+1

    interp_array = np.zeros((len(file_list),nrow, ncol), dtype=np.float32)

    for i,ifile in enumerate(file_list):

        with open(ifile,'r') as f:
            lines = f.readlines()
        f.close()

        # asc file characteristics are gattered from first 6 lines
        # either dx/dy or cellsize should be present
        # NODATA_value may be part of the file specs

        g_nrow = 0
        g_ncol = 0
        g_xllcorner = 0
        g_yllcorner = 0
        g_cellsize = 0
        g_dx = 0
        g_dy = 0
        g_nodata_value = -999999

        for j in range(6):
            line_items = lines[j].split()
            if line_items[0].lower()=='nrows':
                g_nrow = int(line_items[1])
            if line_items[0].lower()=='ncols':
                g_ncol = int(line_items[1])
            if line_items[0].lower()=='xllcorner':
                g_xllcorner = float(line_items[1])
            if line_items[0].lower()=='yllcorner':
                g_yllcorner = float(line_items[1])
            if line_items[0].lower()=='cellsize':
                g_cellsize = float(line_items[1])
            if line_items[0].lower()=='nodata_value':
                g_nodata_value = float(line_items[1])
            if line_items[0].lower()=='dx':
                g_dx= float(line_items[1])
            if line_items[0].lower()=='dy':
                g_dy= float(line_items[1])
        
        if skiprows==None:
            skiprows = 6
        
        asc_grid = np.reshape(np.loadtxt(ifile, skiprows=skiprows),g_nrow*g_ncol)

        g_xy = np.empty((g_nrow*g_ncol, 2), dtype=float)

        inode=0

        if g_cellsize > 0:
            g_dx = g_cellsize
            g_dy = g_cellsize

        for irow in range(g_nrow):
            for icol in range(g_ncol):
                g_xy[inode,0]=g_xllcorner+g_dx*(icol+1/2)
                g_xy[inode,1]=g_yllcorner+g_nrow*g_dy-g_dy*(irow+1/2)
                inode=inode+1

        inter_data = interpolate.griddata(g_xy,asc_grid,model_xy,'nearest')

        inode=0
        for irow in range(nrow):
            for icol in range(ncol):
                if i>0:
                    if inter_data[inode] <=g_nodata_value:
                        interp_array[i,irow,icol]=interp_array[i-1,irow,icol]-1
                    else:
                        interp_array[i,irow,icol]=inter_data[inode]
                else:
                    interp_array[i,irow,icol]=inter_data[inode]
                inode=inode+1

    return interp_array

def Constru_SP(year_0, month_0, day_0, year_F, month_F, day_F, steady_0, steady, timestep, scale='m', n_scale=1, factor=2):
    '''stress period output generator from specific dates and scales

    Args:
        year_0 ('int'): year associated with start date
        month_0 ('int'): month associated with start date
        day_0 ('int'): day associated with start date
        year_F ('int'): year associated with the end date
        month_F ('int'): month associated with the end date
        day_F ('int'): day associated with the end date
        scale ('str'): time window scale e.g. 'd' for scale in days, 'm' for scale in months and 'y' for scale in years (default is 'm')
        n_scale ('int'): step used in the scale (default is 1)
        steady_0 ('bool'): true or False indicating whether or not the first stress period is steady state (default is False)
        timestep ('int'): division of each stress period (default is 18)
        factor ('int'): timestep growth parameter (default is 2)

    Return:
        perlen ('list'): stress period lengths
        fecha_inicioSP ('list'): start date associated with each stress period
        fecha_finSP ('list'): end date associated with each stress period
        nstp ('list'): number of time steps in each stress period
        tsmult ('list'): time step multiplier
        steady ('list'): true or false indicating whether or not stress period is steady state
    '''
        
    # Define empty lists
    perlen = []
    nstp = []
    tsmult = []
    steady = []
    fecha_inicioSP = []
    fecha_finSP = []
    steady = steady if steady is not None else []
    
    # Start and end date with the auxiliary variables
    fecha_0 = datetime.date(year_0,month_0,day_0)
    fecha_F = datetime.date(year_F,month_F,day_F)
    fecha_aux1 = fecha_0
    fecha_aux2 = fecha_0
    
    # Check possible errors with the dates
    if (year_0 or year_F) > 10000 and (year_0 or year_F) > 1000:
        return 'Años son mayores a 10000 o menores a 1000'
        pass
    
    if (fecha_0) >= (fecha_F):
        return 'Fecha inicio mayor que fecha final'
        pass
    
    # The list is made on the scale used
    if scale == 'y': # Yearly scale
        while fecha_aux2<fecha_F:
            fecha_aux2=datetime.date(getattr(fecha_aux2,'year')+n_scale,getattr(fecha_aux2,'month'),getattr(fecha_aux2,'day'))
            if fecha_aux2 > fecha_F:
                fecha_aux2 = fecha_F            
            perlen.append((fecha_aux2-fecha_aux1)/datetime.timedelta(days=1))
            fecha_inicioSP.append(fecha_aux1.isoformat())
            fecha_finSP.append(fecha_aux2.isoformat())
            fecha_aux1 = fecha_aux2
            nstp.append(int(timestep))
            tsmult.append(float(factor))
            steady.append(steady[i] if i < len(steady) else False)
                
    elif scale == 'm': # Monthly scale
        while fecha_aux2 < fecha_F:
            if getattr(fecha_aux2,'month')+n_scale in [13,14,15,16,17,18,19,20,21,22,23]:
                fecha_aux2 = datetime.date(getattr(fecha_aux2,'year')+1,getattr(fecha_aux2,'month')+n_scale-12,getattr(fecha_aux2,'day'))
                if fecha_aux2 > fecha_F:
                    fecha_aux2 = fecha_F                
                perlen.append((fecha_aux2-fecha_aux1)/datetime.timedelta(days=1))
                fecha_inicioSP.append(fecha_aux1.isoformat())
                fecha_finSP.append(fecha_aux2.isoformat())
                fecha_aux1 = fecha_aux2
                nstp.append(int(timestep))
                tsmult.append(float(factor))
                steady.append(steady[i] if i < len(steady) else False)
                
            else:
                fecha_aux2 = datetime.date(getattr(fecha_aux2,'year'),getattr(fecha_aux2,'month')+n_scale,getattr(fecha_aux2,'day'))
                if fecha_aux2 > fecha_F:
                    fecha_aux2 = fecha_F                 
                perlen.append((fecha_aux2-fecha_aux1)/datetime.timedelta(days=1))
                fecha_inicioSP.append(fecha_aux1.isoformat())
                fecha_finSP.append(fecha_aux2.isoformat())
                fecha_aux1 = fecha_aux2
                nstp.append(int(timestep))
                tsmult.append(float(factor))
                steady.append(steady[i] if i < len(steady) else False)

    elif scale == 'd': # Daily scale
        while fecha_aux2 < fecha_F:
            fecha_aux2 += datetime.timedelta(days=n_scale)
            if fecha_aux2 > fecha_F:
                fecha_aux2 = fecha_F
            perlen.append((fecha_aux2-fecha_aux1)/datetime.timedelta(days=1))
            fecha_inicioSP.append(fecha_aux1.isoformat())
            fecha_finSP.append(fecha_aux2.isoformat())
            fecha_aux1 = fecha_aux2
            nstp.append(int(timestep))
            tsmult.append(float(factor))
            steady.append(False)    

    else:
        pass

    # Add if the first stress period is steady
    steady[0] = steady_0
    
    # Stress period outputs
    return perlen, fecha_inicioSP, fecha_finSP, nstp, tsmult, steady

def CI_transporte(sim,CI_block):

    # Solo algunas quimicas
    blq = CI_block[['X','Y','Z','li_int_c','so4_int_','ca_int_c','k_int_co']]
    # Eliminar valores -99.0
    blq = blq[(blq['li_int_c'] != -99.0) & (blq['so4_int_'] != -99.0) & (blq['ca_int_c'] != -99.0) & (blq['k_int_co'] != -99.0)].reset_index(drop=True)

    # Grilla
    model_x = sim.modelgrid.xcellcenters
    model_y = sim.modelgrid.ycellcenters
    model_z = sim.modelgrid.zcellcenters

    nrow = sim.nrow
    ncol = sim.ncol
    nlay = sim.nlay

    # CI a rellenar
    Li_eval = np.zeros((nlay,nrow,ncol))
    SO4_eval = np.zeros((nlay,nrow,ncol))
    Ca_eval = np.zeros((nlay,nrow,ncol))
    K_eval = np.zeros((nlay,nrow,ncol))

    # Coordinadas y quimicas modelo bloque
    coord = blq[['X','Y','Z']].to_numpy()
    Li_blq = blq['li_int_c'].to_numpy()
    SO4_blq = blq['so4_int_'].to_numpy()
    Ca_blq = blq['ca_int_c'].to_numpy()
    K_blq = blq['k_int_co'].to_numpy()

    # Distancia nodos
    kdtree = cKDTree(coord)

    for irow in range(nrow):
        for icol in range(ncol):
            for ilay in range(nlay):
                p0 = (model_x[irow,icol],model_y[irow,icol],model_z[ilay,irow,icol]) 
                dist,index=kdtree.query(p0,1) 
                Li_eval[ilay,irow,icol] = Li_blq[index]
                SO4_eval[ilay,irow,icol] = SO4_blq[index]
                Ca_eval[ilay,irow,icol] = Ca_blq[index]
                K_eval[ilay,irow,icol] = K_blq[index]

    Li_eval = Li_eval.reshape(nlay*nrow*ncol)
    SO4_eval = SO4_eval.reshape(nlay*nrow*ncol)
    Ca_eval = Ca_eval.reshape(nlay*nrow*ncol)
    K_eval = K_eval.reshape(nlay*nrow*ncol)

    return Li_eval,SO4_eval,Ca_eval,K_eval
    

def geology_to_model(sim,geol_units_df):
    '''

    Args:
        sim ('model object'): The model object of type flopy.modflow.mf.Modflow
        geol_units_df ('DataFrame'): Geology dataframe with X,Y,Z and geology index

    Return:
        geol_eval1('array'): Numpy array with the geology evaluated in the model
    '''

    #Coordenadas modelo
    model_x = sim.modelgrid.xcellcenters
    model_y = sim.modelgrid.ycellcenters
    model_z = sim.modelgrid.zcellcenters

    nrow = sim.nrow
    ncol = sim.ncol
    nlay = sim.nlay

    geol_eval = np.zeros((nlay,nrow,ncol),dtype=int)

    coord = geol_units_df[['X','Y','Z']].to_numpy()
    lbl = geol_units_df['GID1'].to_numpy()
    kdtree = cKDTree(coord)

    for irow in range(nrow):
        for icol in range(ncol):
            for ilay in range(nlay):
                p0 = (model_x[irow,icol],model_y[irow,icol],model_z[ilay,irow,icol]) 
                dist,index=kdtree.query(p0,1) 
                geol_eval[ilay,irow,icol] = lbl[index]
    return geol_eval



def check_layer_elevs(elev_array, threshold=None):
    """Check overlaping elevations between layers from a 3D array (nlay,nrow,ncol)

    Args:
        elev_array (3D array): 3D array containing the layer elevation information
        threshold (float, optional): elevation difference assigned if an overlapping layer is encountered. Default value is 0.1.

        Return: 3D array
    """

    if threshold == None:
        threshold = 0.1

    nlay = elev_array.shape[0]
    nrow = elev_array.shape[1]
    ncol = elev_array.shape[2]

    for i in range(1, nlay):
        for irow in range(nrow):
            for icol in range(ncol):
                diff = elev_array[i-1, irow, icol] - elev_array[i, irow, icol]
                if diff < threshold:
                    elev_array[i, irow, icol] = elev_array[i-1, irow, icol] - threshold - 0.00001  # Add a small buffer

    return elev_array

def fill_botm_array(botm,elev_array,gls_lays):
    """fill model botm structured array based on an elev_array

    Args:
        botm (3D array): array to be filled
        elev_array (3D array): elevation array to use for interpolation
        gls_lays (2D array): array containing the number of layers between

    Return:
        botm('array'): 3D array
    """
    nrow = botm.shape[1]
    ncol = botm.shape[2]

    if len(gls_lays) != len(elev_array) - 1:
        raise ValueError("gls_lays should have a length of len(elev_array) - 1")
    
    for irow in range(nrow):
        for icol in range(ncol):
            ibot = elev_array[0,irow,icol]
            j=0
            for i in range(1,len(elev_array)):
                delta = (elev_array[i-1,irow,icol]-elev_array[i,irow,icol])/(gls_lays[i-1])
                for ilay in range(gls_lays[i-1]):
                    ibot = ibot-delta
                    botm[j,irow,icol]=ibot
                    j=j+1

    return botm

def fill_topm_array(botm,interp_elevs,gls_lays):
    """fill model topm structured array based on a botm array

    """

    nlay = botm.shape[0]
    nrow = botm.shape[1]
    ncol = botm.shape[2]

    topm = np.zeros((nlay, nrow, ncol), dtype=float)
    topm[0,:,:] = interp_elevs[0,:,:]
    for i in range(nlay-1):
        topm[i+1,:,:] = botm[i,:,:]
    
    return topm
    
def active_cell(sim,shp_noflow_obj,perc=0.05):
    '''Array generator with active cells

    Args:
        sim ('model object'): The model object of type flopy.modflow.mf.Modflow
        shp_noflow_obj('shapefile'): Shapefile.Reader with non-active cells

    Return:
        noflow('array'): 3D array with 1s or 0s indicating active cells
    '''    
    nrow = sim.nrow
    ncol = sim.ncol
    nlay = sim.nlay

    delr = sim.dis.delr.get_value()
    delc = sim.dis.delc.get_value()

    # Se inicializa la GridIntersect
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)
    
    # Celdas intersectadas y areas
    shapeRecs = shp_noflow_obj.shapeRecords()[0]
    shape = shapeRecs.shape
    p = Polygon(shell=shape.points)
    cells = ix.intersect(p)['cellids']
    areas = ix.intersect(p)['areas']

    # Rellenar noflow matrix
    noflow = np.zeros((nrow,ncol), dtype=int) 
    for id,icell in enumerate(cells):
        if areas[id] < perc*delr*delc:
            pass
        else:
            noflow[icell[0],icell[1]]=1  #reemplazar 0 por ilayer si cada layer tiene 1 shape
        
    # Replicar para todas las capas
    noflow = np.repeat(noflow[None,...], nlay, axis=0)

    return noflow


def strt_head(sim,file_list,ID_UH_array,str_meth='points'):

    if str_meth in ['points','raster']:

        nrow = sim.nrow
        ncol = sim.ncol
        nlay = sim.nlay

        interdata = np.ones((len(file_list),nrow,ncol))*2300.0

        # Se genera un arreglo con las coordenadas de las celdas
        model_xy = np.empty((nrow*ncol, 2), dtype=float)
        model_xy = (sim.modelgrid.xcellcenters,sim.modelgrid.ycellcenters)

        #Se recorre cada archivo para generar las interpolaciones UA-UB-UC a la grilla
        for ifile, file in enumerate(file_list):

            if str_meth == 'points':
                x_test=pd.read_csv(file).loc[:,'x'].to_numpy()
                y_test=pd.read_csv(file).loc[:,'y'].to_numpy()
                values=pd.read_csv(file).loc[:,'VALUE'].to_numpy()
                interdata[ifile] = interpolate.griddata((x_test,y_test),values,model_xy,'nearest').reshape(nrow,ncol)
            if str_meth == 'raster':
                starh_UH = Raster.load(file)
                dem_data = starh_UH.resample_to_grid(sim.modelgrid, band=starh_UH.bands[0], method="nearest", extrapolate_edges=True)
                interdata[ifile] = np.where(dem_data == np.min(dem_data),2300.0,dem_data) #dem_data to numpy, replace -99999.0 

        #Se genera matriz de starting head (nlay,nrow,ncol)
        strt = np.ones(ID_UH_array.shape)*2300.0
        #Se asigna el starting head UA-UB-UC segun el id_UH
        for k in range(nlay):
            strt[k][ID_UH_array[k,:,:] == 0]=interdata[0][ID_UH_array[k,:,:] == 0]
            strt[k][ID_UH_array[k,:,:] == 1]=interdata[1][ID_UH_array[k,:,:] == 1]
            strt[k][ID_UH_array[k,:,:] == 2]=interdata[2][ID_UH_array[k,:,:] == 2]

        return strt

def bc_lat_rch(sim,shp_reclat_obj,q_rl,wel_dtype,criterio_K=0.01,perc=0.15):

    shp_reclat_rec = shp_reclat_obj.records()
    shp_reclat_shapes = shp_reclat_obj.shapes()
    noflow=sim.bas6.ibound.array
    shapeRecs = shp_reclat_obj.shapeRecords()   
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)

    #Se inician variables que guardaran los nodos de RecLat, caudal, id recarga, Layer
    nodos_reclat = []
    qin_reclat = []
    rec_reclat = []
    lay_reclat = []
    Li_reclat = []
    SO4_reclat = []
    Ca_reclat = []
    K_reclat = []

    #se recorre cada shape asociado a la polilínea de recarga lateral
    for shaperec in shapeRecs:
        celdas_reclat = []
        nrec = 0

        #Se guarda el id de recarga y layer para cada polilínea con la recarga lateral
        rec = shaperec.record['id']
        from_layer = shaperec.record['from_layer']
        to_layer = shaperec.record['to_layer']

        #Se intersecta la polilínea con la grilla, guardando los id de la celda y el largo de intersección
        pline = shaperec.shape.points
        ls = LineString(pline)
        cells = ix.intersect(ls)['cellids']
        largo_inter = ix.intersect(ls)['lengths']

        delr = sim.dis.delr.get_value()

        if cells.any():
            # Se recorre cada celda intersectada por la polilínea, según los layer en los cuales está la CB
            #para guardar cada una de ellas
            for k in range(from_layer,to_layer+1):
                for id,icell in enumerate(cells):
                    # No se consideran celdas que se encuentran inactivas
                    if noflow[k-1,icell[0],icell[1]] == 0:
                        pass
                    # No se consideran celdas que la intersección sea menor al 'perc' % del largo de la celda
                    elif largo_inter[id] < perc*delr:
                        pass
                    elif sim.lpf.hk.array[k-1][icell[0]][icell[1]] < criterio_K:
                        pass
                    else:
                        #Se guarda el id de recarga, la celda, y el layer de la celda intersectada
                        rec_reclat.append(rec)
                        celdas_reclat.append(icell)
                        lay_reclat.append(k)

            nrec = len(celdas_reclat)

            if nrec>0:
                # se divide el caudal total de la polilínea por la cantidad de celdas intersectadas
                iq_rl = q_rl[q_rl['id_recarga'] == rec]['Q(m3/d)'].tolist()
                iq_rl=iq_rl[0]/nrec
                Li = q_rl[q_rl['id_recarga'] == rec]['Li'].values
                SO4 = q_rl[q_rl['id_recarga'] == rec]['SO4'].values
                Ca = q_rl[q_rl['id_recarga'] == rec]['Ca'].values
                K = q_rl[q_rl['id_recarga'] == rec]['K'].values

                #Se asocia un caudal a cada celda intersectada
                for icell in celdas_reclat:
                    nodos_reclat.append(icell)
                    qin_reclat.append(iq_rl)
                    Li_reclat.append(Li)
                    SO4_reclat.append(SO4)
                    Ca_reclat.append(Ca)
                    K_reclat.append(K)

                nrec = len(celdas_reclat)

    sp1 = np.zeros((len(nodos_reclat)), dtype=wel_dtype)
    sp1 = sp1.view(np.recarray)

    # se leen las celdas, y se asigna una variable auxiliar asociada con la zona de recarga
    if len(wel_dtype) == 5:
        for icell in range(len(nodos_reclat)):
            sp1[icell] = (lay_reclat[icell]-1,nodos_reclat[icell][0],nodos_reclat[icell][1]
                      ,qin_reclat[icell],rec_reclat[icell]+100)
    else:
        for icell in range(len(nodos_reclat)):
            sp1[icell] = (lay_reclat[icell]-1,nodos_reclat[icell][0],nodos_reclat[icell][1]
                      ,qin_reclat[icell],rec_reclat[icell]+100,Li_reclat[icell],SO4_reclat[icell],Ca_reclat[icell],K_reclat[icell])

    return sp1


def well_to_day(df):
    '''Convierte del DataFrame mensual en un caudal uniforme diario
    retornando un diccionario con los caudales diarios de todos los pozos'''
    
    dic_well = {} #diccionario vacio donde se almacenará el DataFrame de cada pozo
    name_well = df['pozo'].unique() #nombre de los pozos
    num_well = len(df['pozo'].unique()) #cantidad de pozos

    for i in range(num_well): #recorrer los pozos
        df_well = df[df['pozo'] == name_well[i]]
        df_well.reset_index(drop=True,inplace=True)  

        df_new = pd.DataFrame(columns = ['Days','Flow']) #DataFrame a rellenar

        for j in range(len(df_well['Fecha_inicio'])): #recorrer meses de cada pozo
            date_i = datetime.datetime.strptime(df_well['Fecha_inicio'][j], "%d/%m/%Y")
            date_f = datetime.datetime.strptime(df_well['Fecha_fin'][j], "%d/%m/%Y")
            delta = date_f - date_i

            days = []
            for k in range(delta.days): #agregar dias del mes
                days.append(date_i + datetime.timedelta(days=k))

            #almecenar datos
            df_aux = pd.DataFrame()
            df_aux['Days'] = days
            df_aux['Flow'] = np.repeat(df_well['Caudal_m3/d'][j],len(days))

            #concaquetar datos
            df_new = pd.concat([df_new,df_aux],axis=0)

        #guardar datos en diccionario
        df_new.reset_index(drop=True,inplace=True)
        dic_well[name_well[i]] = df_new
    
    df_q = pd.DataFrame(columns = ['Borehole','Days','Flow']) #DataFrame a rellenar
    for key, value in dic_well.items():
        value['Borehole'] = np.repeat(key,len(value))
        df_q = pd.concat([df_q,value],axis=0)
    df_q.reset_index(drop=True,inplace=True)

    return df_q

def well_to_SP(SP,well_q):
    
    #Data del SP
    fecha_inicialSP = SP[1]
    fecha_finalSP = SP[2]
        
    df_sp = pd.DataFrame(columns = ['Borehole','sp','Fecha_inicio','Fecha_fin','Caudal']) #DataFrame a rellenar

    for group, fname in well_q.groupby('Borehole'):
        sp = []
        Caudal = []
        for i in range(len(fecha_inicialSP)):
            sp.append(i)

            #Filtrar segun fechas
            f_i = fname['Days'] >= fecha_inicialSP[i]
            f_f = fname['Days'] < fecha_finalSP[i]
            filt = f_i & f_f
            Caudal.append(np.mean(fname.loc[filt]['Flow'])) #Caudal promedio segun fechas establecidas

        dic = {'Borehole':[group]*len(sp),'sp':sp,'Fecha_inicio':fecha_inicialSP,'Fecha_fin':fecha_finalSP,'Caudal':Caudal}         
        df_sp_well = pd.DataFrame(dic) #DataFrame de un pozo
            
        df_sp = pd.concat([df_sp,df_sp_well],axis=0) 
            

    df_sp.dropna(subset = ["Caudal"], inplace=True) #Eliminar filas NaN
    df_sp.reset_index(drop=True,inplace=True)       #Reiniciar index

    # Agregar 0 para el ultimo stress periods si la data es menor a todos los SP
    df_fname = pd.DataFrame(columns = df_sp.columns)
    for group, fname in df_sp.groupby('Borehole'): 
        if fname.iloc[-1]['sp']<len(SP[1])-1:
            fname.loc[len(fname.index),df_sp.columns] = fname.iloc[-1]
            fname.loc[len(fname.index)-1,'sp'] = fname.iloc[-1]['sp'] +1
            fname.loc[len(fname.index)-1,'Fecha_inicio'] = (pd.to_datetime(fname.iloc[-1]['Fecha_inicio']) + pd.DateOffset(months=1)).strftime('%Y-%m-%d')
            fname.loc[len(fname.index)-1,'Fecha_fin'] = (pd.to_datetime(fname.iloc[-1]['Fecha_fin']) + pd.DateOffset(months=1)).strftime('%Y-%m-%d')
            fname.loc[len(fname.index)-1,'Caudal'] = 0
            df_fname = pd.concat([df_fname,fname])
        else: 
            df_fname = pd.concat([df_fname,fname])
    df_fname.reset_index(drop=True,inplace=True)

    return df_fname

def assign_layer(sim,wells_ubi, iny = False):

    z_bot = sim.modelgrid.botm #bot de toda la grilla
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True) #intersección dimensiones y celdas

    # Filtrar aquellos que no intersectan
    wells_ubi['intersecta']=True
    for index,row in wells_ubi.iterrows():
        try:
            ix.intersects(Point(row['UTM_East'],row['UTM_North']))['cellids'][0][0]
        except:
            wells_ubi.loc[index,'intersecta']=False
    wells_ubi=wells_ubi[wells_ubi['intersecta']==True]
    wells_ubi.reset_index(inplace=True,drop=True)

    
    if len(wells_ubi)>0:
        wells_ubi['CELL_X'] = wells_ubi.apply(lambda x: ix.intersects(Point(x['UTM_East'],x['UTM_North']))['cellids'][0][0], axis=1) #row
        wells_ubi['CELL_Y'] = wells_ubi.apply(lambda x: ix.intersects(Point(x['UTM_East'],x['UTM_North']))['cellids'][0][1], axis=1) #col

        if iny:
            #asignar los layer a cada uno de los pozos, siendo igual a 1 (son pozos someros)
            wells_ubi['top_layer'] = 0          
            wells_ubi['bot_layer'] = 0
        
        else: 
            wells_ubi['eval_depth'] = wells_ubi.apply(lambda x: z_bot[:,x['CELL_X'],x['CELL_Y']], axis=1) #bot de todas las celdas para un row col seleccionado
            # seleccionar el top de la celda mas cercana dependiendo del signo 
            wells_ubi['dif_top'] = wells_ubi.apply(lambda x: x['Top_scr'] - x['eval_depth'],axis=1)     
            wells_ubi['min_top'] = wells_ubi.apply(lambda x: np.argmin(np.abs(x['dif_top'])),axis=1)   
            wells_ubi['sign_top'] = wells_ubi.apply(lambda x: np.sign(x['dif_top'][x['min_top']]),axis=1)   
            wells_ubi['top_layer'] = wells_ubi.apply(lambda x: x['min_top'] if x['sign_top'] > 0 else x['min_top']+1, axis=1) 
            # seleccionar el bot de la celda mas cercana dependiendo del signo 
            wells_ubi['dif_bot'] = wells_ubi.apply(lambda x: x['Bot_scr'] - x['eval_depth'],axis=1)     
            wells_ubi['min_bot'] = wells_ubi.apply(lambda x: np.argmin(np.abs(x['dif_bot'])),axis=1)   
            wells_ubi['sign_bot'] = wells_ubi.apply(lambda x: np.sign(x['dif_bot'][x['min_bot']]),axis=1)   
            wells_ubi['bot_layer'] = wells_ubi.apply(lambda x: x['min_bot'] if x['sign_bot'] > 0 else x['min_bot']+1, axis=1) 

        return wells_ubi

    else:
        return wells_ubi

def dic_sp_well(df_well_to_sp,wells_ubi,wel_dtype):
    stress_period_data = {}
    pozos_iny = ['IN-01','IN-01R','IN-01R2','IN-02','IN-02R','IN-02R2','IN-03','IN-03R','IN-03R2','IN-03R3','IN-04','IN-04R','IN-04R3','IN-05','IN-05R','IN-06','IN-06R','IN-07','IN-07R']

    if len(wells_ubi)>0:

        df_well_to_sp['auxvar'] = df_well_to_sp.apply(lambda x: 50 if x['Borehole'] in pozos_iny else 55, axis=1)
        for group, fname in df_well_to_sp.groupby('sp'):
            #Solo seleccionar aquellos pozos que esten en el df de wells_ubi
            fname2=fname.loc[fname.Borehole.isin(wells_ubi['Borehole'])]
            if len(wel_dtype) == 5: #Sin quimica
                sp = fname2.merge(wells_ubi,how='left')[['bot_layer','CELL_X','CELL_Y','Caudal','auxvar']].to_records(index=False).astype(wel_dtype)
            else: #Con quimica
                sp = fname2.merge(wells_ubi,how='left')[['bot_layer','CELL_X','CELL_Y','Caudal','auxvar','Li_pct','SO4_pct','Ca_pct','K_pct']].to_records(index=False).astype(wel_dtype)
            stress_period_data.update({group:sp})


        return stress_period_data

    else:
        sp1 = np.zeros((len(stress_period_data)), dtype=wel_dtype)
        sp1 = sp1.view(np.recarray)
        stress_period_data={0:sp1}
        
        return stress_period_data

def edit_sms(sim):

    modelname = sim.name
    model_ws = sim.model_ws

    with open(os.path.join(model_ws,modelname+".sms"),'r') as f:
        lines = f.readlines()
    f.close()

    #DAMPBOT no es siempre necesario. No se si se deba dejar esto fijo.
    lines[1]=lines[1][:-1]+' SOLVEACTIVE DAMPBOT\n'
    lines[3]=lines[3][:-5]
    #lines[3]='1 2 7 14 0 0.0 1 0.001 \n'

    with open(os.path.join(model_ws,modelname+".sms"),'w') as g:
        g.writelines(lines)
    g.close()

def edit_oc(sim,compact_budget=False):

    modelname = sim.name
    model_ws = sim.model_ws
    nper = sim.nper

    with open(os.path.join(model_ws,modelname+".oc"),"w") as f:
        f.write(' ATSA NPSTPS 1000\n')
        f.write('  HEAD SAVE UNIT 51\n')
        f.write('  HEAD PRINT FORMAT 0\n')
        #f.write('  DRAWDOWN SAVE UNIT 31\n')
        #f.write('  DRAWDOWN PRINT FORMAT 0\n')
        f.write('  CONC SAVE UNIT 132\n')
        f.write('  CONC PRINT FORMAT 0 \n')
        if compact_budget:
            f.write('  COMPACT BUDGET\n')
        
        for iper in range(nper):
            f.write('PERIOD '+str(iper+1)+'\n')
            #f.write('   DELTAT 1.000000e+00\n')
            #f.write('   TMINAT 1.000000e-01\n')
            #f.write('   TMAXAT 2.000000e+02\n')
            #f.write('   TADJAT 2.000000e+00\n')
            #f.write('   TCUTAT 2.000000e+00\n')
            f.write('   SAVE HEAD\n')
            #f.write('   SAVE DRAWDOWN\n')
            f.write('   SAVE BUDGET\n')
            f.write('   PRINT BUDGET\n')
            #f.write('  SAVE CONC\n')

def edit_lpf(sim,vani_opt=False):

    modelname = sim.name
    model_ws = sim.model_ws
    nlay = sim.nlay

    with open(os.path.join(model_ws,modelname+".lpf"),'r') as f:
        lines = f.readlines()
    f.close()
    with open(os.path.join(model_ws,modelname+".lpf"),'w') as g:
        for line in range(len(lines)):
            g.write(lines[line])
            if line == 6:
                break
        for i in range(nlay):
            g.write("OPEN/CLOSE "+"kx_"+str(i+1)+".dat 1.0 (FREE) -1\n")
            if vani_opt:
                g.write("OPEN/CLOSE "+"vani_"+str(i+1)+".dat 1.0 (FREE) -1\n")
            else:
                g.write("OPEN/CLOSE "+"kz_"+str(i+1)+".dat 1.0 (FREE) -1\n")
            g.write("OPEN/CLOSE "+"ss_"+str(i+1)+".dat 1.0 (FREE) -1\n")
            g.write("OPEN/CLOSE "+"sy_"+str(i+1)+".dat 1.0 (FREE) -1\n")
    f.close()

def areas_evt(sim,shp_rdir_obj,zone_rch_infiltracion):

    noflow=sim.bas6.ibound.array[0]
    shapeRecs = shp_rdir_obj.shapeRecords()
    nrow = sim.nrow
    ncol = sim.ncol
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)

    #arreglo que contendrá que cantidad de cada zona de EVT (en área) hay en cada celda
    areas_arr = np.zeros((len(shapeRecs),nrow,ncol), dtype=float) 
    areas_arr_id = np.zeros((len(shapeRecs),nrow,ncol), dtype=float)

    for ishape,shaperecs in enumerate(shapeRecs):
        id_EVT=shaperecs.record['id']
        shape = shaperecs.shape 
        polyg_points = shape.points #pasar a puntos
        p = Polygon(shell=polyg_points) #volver a polígonos
        result = ix.intersect(p) #intersectar 
        cellids = result['cellids']
        areas = result['areas']
        for icell in range(len(cellids)):
            if noflow[cellids[icell][0],cellids[icell][1]] == 0:
                areas_arr[ishape,cellids[icell][0],cellids[icell][1]] = 0
                areas_arr_id[ishape,cellids[icell][0],cellids[icell][1]]=0
            else:
                areas_arr[ishape,cellids[icell][0],cellids[icell][1]] = areas[icell]
                areas_arr_id[ishape,cellids[icell][0],cellids[icell][1]]=id_EVT

    #arreglo que contendrá la zona de EVT final en cada celda 
    areas_arr_s = np.zeros((nrow,ncol), dtype=int) 

    #arreglo que contendrá la cantidad de area de las zonas de EVT para cada celda 
    areas = np.zeros(areas_arr.shape[0], dtype=int) 

    #se asigna la correspondencia entre las zonas de los SHP (ipol), con las zonas de EVT (mayor área en la celda)
    for irow in range(nrow):
        for icol in range(ncol):
            for ipol in range(areas_arr.shape[0]):
                areas[ipol] = areas_arr[ipol,irow,icol]
            areas_arr_s[irow,icol] = areas_arr_id[np.argmax(areas),irow,icol]*noflow[irow,icol]

    # se asigna zona 1 a aquellas que tienen infiltración de pozas
    areas_arr_s[zone_rch_infiltracion==1]=1 
                
    return areas_arr_s



def write_ETS(sim,ETS_series,areas_arr_s,seg=6):

    modelname = sim.name
    model_ws = sim.model_ws
    nrow = sim.nrow
    ncol = sim.ncol
    nper = sim.nper
    top = sim.dis.gettop() #top layer1

    #Escribir paquete ETS. forma extendida
    with open(os.path.join(model_ws,modelname+".ets"),"w") as f:
        f.write("# MODFLOW-USGs Evapotranspiration Package\n")
        f.write(" 1 50 0  6  0\n")

        cont=1
        for i in range(nper):
            f.write("  1 1 1 -1 1\n")
            f.write("INTERNAL  1.000000e+000  (FREE)  -1  ET Surface\n")
            cont=1
            for i in range(nrow):
                for j in range(ncol-1):
                    a="{:e}".format(top[i][j])
                    f.write(str(a)+" ")
                    if (cont)%10==0:
                        f.write("\n")
                    cont=cont+1
                a="{:e}".format(top[i][j])
                f.write(str(a)+" ")                    
                
                f.write("\n")
                cont=1
            #f.write("\n")
        
            f.write("INTERNAL  1.000000e+000  (FREE)  -1  ET Rate\n")
            cont=1
            param='ET Rate'
            df_aux=ETS_series.loc[ETS_series.Parameter==param,:].drop('Parameter',axis=1).to_numpy()
            for i in range(nrow):
                for j in range(ncol-1):
                    #zone='Zone'+str(areas_arr_s[i,j])
                    a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                    f.write(str(a)+" ")
                    if (cont)%10==0:
                        f.write("\n")
                    cont=cont+1
                a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                f.write(str(a)+" ")

                f.write("\n")
                cont=1
            #f.write("\n")

            f.write("INTERNAL  1.000000e+000  (FREE)  -1  Extinction Depth\n")
            cont=1
            param='Extinction Depth'
            df_aux=ETS_series.loc[ETS_series.Parameter==param,:].drop('Parameter',axis=1).to_numpy()
            for i in range(nrow):
                for j in range(ncol-1):
                    a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                    f.write(str(a)+" ")
                    if (cont)%10==0:
                        f.write("\n")
                    cont=cont+1
                a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                f.write(str(a)+" ")
                
                f.write("\n")
                cont=1
            #f.write("\n")
            for iseg in range(seg-1):
                f.write("INTERNAL  1.000000e+000  (FREE)  -1  PXDP"+str(iseg)+"\n")
                cont=1
                param='PXDP'+str(iseg)
                df_aux=ETS_series.loc[ETS_series.Parameter==param,:].drop('Parameter',axis=1).to_numpy()
                for i in range(nrow):
                    for j in range(ncol-1):
                        a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                        f.write(str(a)+" ")
                        if (cont)%10==0:
                            f.write("\n")
                        cont=cont+1
                    a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                    f.write(str(a)+" ")
                    
                    f.write("\n")
                    cont=1
                #f.write("\n")

                f.write("INTERNAL  1.000000e+000  (FREE)  -1  PETM"+str(iseg)+"\n")
                cont=1
                param='PETM'+str(iseg)
                df_aux=ETS_series.loc[ETS_series.Parameter==param,:].drop('Parameter',axis=1).to_numpy()
                for i in range(nrow):
                    for j in range(ncol-1):
                        a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                        f.write(str(a)+" ")
                        if (cont)%10==0:
                            f.write("\n")
                        cont=cont+1
                    a="{:e}".format(df_aux[0,areas_arr_s[i,j]])
                    f.write(str(a)+" ")
            
                    f.write("\n")
                    cont=1
                #f.write("\n")


    f.close()

def write_rch(sim,tasa_pond,dic_qca=None):
    modelname = sim.name
    model_ws = sim.model_ws
    nrow = sim.nrow
    ncol = sim.ncol
    nper = sim.nper

    if dic_qca == None: #Sin quimica
        with open(os.path.join(model_ws,modelname+".rch"),"w") as f:
            f.write("# MODFLOW-USGs Recharge Package\n")
            f.write(" 3 50\n")
            for per in range(nper):
                f.write(" 1\n") #si tiene concentracion, hay que poner CONC al lado de este 1, y desúes de la tasa de recarga, agregar la concentración con el mismo formato
                f.write("INTERNAL  1.000000e+000  (FREE)  -1  RECHARGE\n")
                cont=1
                for i in range(nrow):
                    for j in (range(ncol-1)):
                        a="{:e}".format(tasa_pond[per][i][j])
                        f.write(str(a+" "))
                        if cont % 10==0:
                            f.write("\n")
                        cont=cont+1
                    a="{:e}".format(tasa_pond[per][i][j+1])
                    f.write(str(a+" "))

                    f.write("\n")
                    cont=1
                # f.write("\n")
            f.close()

    else: #Con quimicas
        with open(os.path.join(model_ws,modelname+".rch"),"w") as f:
            f.write("# MODFLOW-USGs Recharge Package\n")
            f.write(" 3 50 CONC\n")
            f.write(" 1 1 1 1\n")
            for per in range(nper):
                f.write(" 1 INCONC\n") #si tiene concentracion, hay que poner CONC al lado de este 1, y desúes de la tasa de recarga, agregar la concentración con el mismo formato
                f.write("INTERNAL  1.000000e+000  (FREE)  -1  RECHARGE\n")
                cont=1
                for i in range(nrow):
                    for j in (range(ncol-1)):
                        a="{:e}".format(tasa_pond[per][i][j])
                        f.write(str(a+" "))
                        if cont % 10==0:
                            f.write("\n")
                        cont=cont+1
                    a="{:e}".format(tasa_pond[per][i][j+1])
                    f.write(str(a+" "))

                    f.write("\n")
                    cont=1
                # f.write("\n")
                # QUIMICAS
                for value in dic_qca.values():
                    f.write("INTERNAL  1.000000e+000  (FREE)  -1  RECHARGE CONC\n")
                    cont=1
                    for i in range(nrow):
                        for j in (range(ncol-1)):
                            a="{:e}".format(value[per][i][j])
                            f.write(str(a+" "))
                            if cont % 10==0:
                                f.write("\n")
                            cont=cont+1
                        a="{:e}".format(value[per][i][j+1])
                        f.write(str(a+" "))

                        f.write("\n")
                        cont=1
                    # f.write("\n")
            f.close()

def zonas_pozas(sim,filename,nfield):
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)
    shp_rdir_obj = shp.Reader(filename)
    shapeRecs = shp_rdir_obj.shapeRecords()

    nrow = sim.nrow
    ncol = sim.ncol
    zone_array = np.zeros((nrow,ncol))

    areas_arr = np.zeros((len(shapeRecs),nrow,ncol), dtype=float) #arreglo que contendrá que cantidad de cada zona de
    a=int(nfield)                                                              #recarga (en área) hay en cada celda

    zonas=np.zeros((1,areas_arr.shape[0]), dtype=float) 
    
    for ishape in range(len(shapeRecs)):
        shape = shapeRecs[ishape].shape
        rec_zone = shapeRecs[ishape].record[:a][a-1] ######## verificar que se puede leer la zona --> si funciona pero lee el 1° campo
        zonas[0,ishape] = int(rec_zone)  #tenemos la zona de cada poligono. verificar si se lee la zona --> si funciona
        polyg_points = shape.points #pasar a puntos
        p = Polygon(shell=polyg_points) #volver a polígonos
        result = ix.intersect(p) #intersectar 
        cellids = result['cellids']
        areas = result['areas']
        for icell in range(len(cellids)):
            areas_arr[ishape,cellids[icell][0],cellids[icell][1]]=areas[icell]

    areas = np.zeros((1,areas_arr.shape[0]), dtype=float) #arreglo q contendrá la cantidad de area de las zonas de recarga 
                                                        #para cada celda 
    for irow in range(nrow):
        for icol in range(ncol):
            for ipol in range(areas_arr.shape[0]):
                areas[0,ipol] = areas_arr[ipol,irow,icol]
            zone_array[irow,icol]=zonas[0,np.argmax(areas,axis=1)]   ###verificar si funciona --> si funciona
            #zona_array es un arreglo que para cada celda le asigna una zona. ojo q el numerativo del array empieza en 0

    return zone_array, areas_arr

def recarga_mensual(filename1,filename2):
    # Calculamos la recarga mensual, la pasamos a m3/d dividiendo por 86.4 y dividimos tb por el area total y días del mes
    # filename1 es un csv con la fecha y la infiltración diaria en l/s
    # filename2 es un csv con el área de cada zona
    df = pd.DataFrame(data=filename1)

    df.Fecha = pd.to_datetime(df.Fecha, format='%d/%m/%Y')
    df.set_index('Fecha', inplace=True)
    
    df_mes = df.resample('MS').sum()  #calculo de recarga mensual
    
    #ahora las conversiones
    #primero obtener los dias de cada mes
    df_mes_R = df_mes.reset_index() 
    df_fecha = df_mes_R.iloc[0:,0]
    df_fecha = df_fecha.apply(lambda t: pd.Period(t, freq='S').days_in_month).tolist() #se pasa a lista
    inv_fecha = [1/i for i in df_fecha] #se calcula el inverso para futuros calculos
    
    #segundo obtener el área de cada zona
    areas_rec = filename2['Area_m2'].tolist()
    a = [i * (1/86.4) for i in areas_rec] #se multiplica por un factor para despues poder pasar la recarga de l/s a m3/d
    
    #calculos finales
    t = df_mes/a #se divide por el area total de la zona
    tasa = t.mul(inv_fecha, axis=0) #se multiplica el anterior el inverso de la cantidad de días del mes
    
    return tasa


def rch_nat_sp(SP,inf_diaria,df_nat):

    perlen = SP[0]
    fecha_inicialSP = SP[1]
    fecha_finalSP = SP[2]

    inf_diaria.Fecha = pd.to_datetime(inf_diaria.Fecha, format='%Y-%m-%d')
    inf_diaria.set_index('Fecha', inplace=True)  

    inf_array = np.zeros((len(perlen),len(inf_diaria.columns)))

    for i in range(len(fecha_finalSP)):
        fecha_final = datetime.datetime.strptime(fecha_finalSP[i], "%Y-%m-%d") - datetime.timedelta(days=1) #ARREGLAR
        inf_array[i,:] = np.sum(inf_diaria[fecha_inicialSP[i]:fecha_final],axis=0)

    df_rch_sp = pd.DataFrame(data = inf_array, index = SP[1], columns = inf_diaria.columns)
    df_rch_sp[inf_diaria.columns] = df_rch_sp.mul((1/np.array(perlen)),axis=0)

    df_rch_sp.reset_index(inplace=True)

    return df_rch_sp

def rch_poza_acop_sp(SP,inf_diaria,df_poza_acop):

    perlen = SP[0]
    fecha_inicialSP = SP[1]
    fecha_finalSP = SP[2]

    inf_diaria.Fecha = pd.to_datetime(inf_diaria.Fecha, format='%Y-%m-%d')
    inf_diaria.set_index('Fecha', inplace=True)  

    inf_array = np.zeros((len(perlen),len(inf_diaria.columns)))

    for i in range(len(fecha_finalSP)):
        fecha_final = datetime.datetime.strptime(fecha_finalSP[i], "%Y-%m-%d") - datetime.timedelta(days=1) #ARREGLAR
        inf_array[i,:] = np.sum(inf_diaria[fecha_inicialSP[i]:fecha_final],axis=0)

    df_rch_sp = pd.DataFrame(data = inf_array, index = SP[1], columns = inf_diaria.columns)
    df_rch_sp[inf_diaria.columns] = df_rch_sp.mul((1/np.array(perlen)),axis=0).mul((1/df_poza_acop.A_Name.to_numpy()),axis=1) 

    df_rch_sp.reset_index(inplace=True)

    return df_rch_sp

def inf_to_day(filename):

    # DF infiltración mensual
    df_inf = pd.DataFrame(data=filename)

    df_new = pd.DataFrame(columns = ['Days','Inf','MOP','SOP']) #DataFrame a rellenar

    for i in range(len(df_inf)-1):
        date_i = datetime.datetime.strptime(df_inf['Meses'][i], "%d/%m/%Y")
        date_f = datetime.datetime.strptime(df_inf['Meses'][i+1], "%d/%m/%Y")
        delta = date_f - date_i

        days = []
        for j in range(delta.days): #agregar dias del mes
            days.append(date_i + datetime.timedelta(days=j))

        df_aux = pd.DataFrame()
        df_aux['Days'] = days
        # Dividir volumen en partes iguales de los días
        df_aux['Inf'] = np.repeat(df_inf['Inf'][i]/len(days),len(days))
        df_aux['MOP'] = np.repeat(df_inf['MOP'][i]/len(days),len(days))
        df_aux['SOP'] = np.repeat(df_inf['SOP'][i]/len(days),len(days))

        # concaquetar datos
        df_new = pd.concat([df_new,df_aux],axis=0)

    df_new.reset_index(drop=True,inplace=True)    

    return df_new

def rch_pozas_sp(SP,filename,i_dic):
    
    #Data del SP
    perlen = SP[0]
    fecha_inicialSP = SP[1]
    fecha_finalSP = SP[2]

    df = pd.DataFrame(data = filename)

    df.Days = pd.to_datetime(df.Days, format='%d/%m/%Y')
    df.set_index('Days', inplace=True) 
    col = df.columns.to_list()
    inf_array = np.zeros((len(perlen),len(col)))

    for i in range(len(fecha_finalSP)):
        fecha_final = datetime.datetime.strptime(fecha_finalSP[i], "%Y-%m-%d") - datetime.timedelta(days=1) #ARREGLAR
        inf_array[i,:] = np.sum(df[fecha_inicialSP[i]:fecha_final],axis=0)

    df_rch_pozas_sp = pd.DataFrame(data = inf_array, index = SP[1], columns = col)
    inv_fecha = [1/i for i in perlen]

    i_count = i_dic['i_count']
    i_mop = i_dic['i_mop']
    i_sop = i_dic['i_sop']

    if i_count==0:
        df_rch_pozas_sp.Inf *=0
    else:
        df_rch_pozas_sp.Inf *=(1/i_count)

    if i_mop==0:
        df_rch_pozas_sp.MOP *=0
    else:
        df_rch_pozas_sp.MOP *=(1/i_mop)

    if i_sop==0:
        df_rch_pozas_sp.SOP *=0
    else:
        df_rch_pozas_sp.SOP *=(1/i_sop)


    tasa_p = df_rch_pozas_sp.mul(inv_fecha, axis=0)

    return tasa_p


def area_array(sim,areas_rch):

    nrow = sim.nrow
    ncol = sim.ncol

    areas_recarga = areas_rch['rch_nat']
    areas_recarga_pozas = areas_rch['rch_pozas']
    areas_recarga_acopio_m = areas_rch['rch_MOP']
    areas_recarga_acopio_s = areas_rch['rch_SOP']

    areas_arr_s_tot = np.zeros((nrow,ncol), dtype=float)
    areas_arr_mop = np.zeros((nrow,ncol), dtype=float)
    areas_arr_sop = np.zeros((nrow,ncol), dtype=float)

    i_count = 0
    i_mop = 0
    i_sop = 0

    for irow in range(nrow):
        for icol in range(ncol):
            iarea = 0
            imop = 0
            isop = 0
            for ipol in range(areas_recarga_pozas.shape[0]):
                iarea += areas_recarga_pozas[ipol,irow,icol]
                i_count += areas_recarga_pozas[ipol,irow,icol]
            for imo in range(areas_recarga_acopio_m.shape[0]):
                imop += areas_recarga_acopio_m[imo,irow,icol]
                i_mop += areas_recarga_acopio_m[imo,irow,icol]
            for iso in range(areas_recarga_acopio_s.shape[0]):
                isop += areas_recarga_acopio_s[iso,irow,icol]
                i_sop += areas_recarga_acopio_s[iso,irow,icol]

            
            areas_arr_s_tot[irow,icol] = iarea
            areas_arr_mop[irow,icol] = imop
            areas_arr_sop[irow,icol] = isop

    areas_dic = {'areas_arr_s_tot':areas_arr_s_tot,
                 'areas_arr_mop':areas_arr_mop,
                 'areas_arr_sop':areas_arr_sop}
    i_dic = {'i_count':i_count,
             'i_mop':i_mop,
             'i_sop':i_sop}

    return areas_dic, i_dic

def tasa_pond_func(sim,zonas_recarga,tasa_nat,tasa_p,areas_dic):

    nrow = sim.nrow
    ncol = sim.ncol
    nper = sim.nper

    delc = sim.dis.delc.get_value()
    delr = sim.dis.delr.get_value()

    Area_nodo = delc * delr

    tasa_pond = np.zeros((nper,nrow,ncol),dtype=float)

    areas_arr_s_tot = areas_dic['areas_arr_s_tot']
    areas_arr_mop = areas_dic['areas_arr_mop']
    areas_arr_sop = areas_dic['areas_arr_sop']
    
    for isp in range (nper):
        for irow in range(nrow):
            for icol in range(ncol):
                zona_rn = int(zonas_recarga[irow,icol])
                itasa = tasa_nat[zona_rn][isp] #tasa recarga natural en el nodo
                iarea_pozas = areas_arr_s_tot[irow,icol] # area total de pozas en el nodo
                iarea_ac_m = areas_arr_mop[irow,icol] # area total de pozas en el nodo
                iarea_ac_s = areas_arr_sop[irow,icol] # area total de pozas en el nodo
                iarea_nat = Area_nodo - iarea_pozas - iarea_ac_m - iarea_ac_s
                itasa_pozas = tasa_p["Inf"][isp] #tasa pozas
                itasa_aco_m = tasa_p["MOP"][isp] #tasa acopios MOP
                itasa_aco_s = tasa_p["SOP"][isp] #tasa acopios MOP            
                itasa_pond = (itasa * iarea_nat + itasa_pozas * iarea_pozas + itasa_aco_m * iarea_ac_m + itasa_aco_s * iarea_ac_s) / Area_nodo
                tasa_pond[isp,irow,icol] = itasa_pond

    return tasa_pond

def tasa_pond_array(sim,df_cell_rch,tasa_nat,tasa_poza,tasa_acopios):

    delc = sim.dis.delc.get_value()
    delr = sim.dis.delr.get_value()
    A_cell = delc*delr

    df_cell_rch['T_Nat'] = df_cell_rch.apply(lambda x: (tasa_nat[list(map(str,x['Nat']))]*(x['A_Nat']/A_cell)).sum(axis=1).to_numpy(),axis=1)
    df_cell_rch['T_Poza'] = df_cell_rch.apply(lambda x: (tasa_poza[x['Poza']]*(x['A_Poza']/A_cell)).sum(axis=1).to_numpy(),axis=1)
    df_cell_rch['T_Acop'] = df_cell_rch.apply(lambda x: (tasa_acopios[x['Acop']]*(x['A_Acop']/A_cell)).sum(axis=1).to_numpy(),axis=1)
    df_cell_rch['T_pond'] = df_cell_rch.apply(lambda x: x['T_Nat'] + x['T_Poza'] + x['T_Acop'],axis=1)

    df_cell_rch.set_index('Cell',inplace=True)

    nrow = sim.nrow
    ncol = sim.ncol
    nper = sim.nper

    tasa_pond = np.zeros((nper,nrow,ncol))

    for i in range(nrow):
        for j in range(ncol):
            tasa_pond[:,i,j] = df_cell_rch.xs((i,j)).T_pond
    
    df_cell_rch.reset_index(inplace=True)

    return tasa_pond

def hydprop_to_model(df_hydprop,geol_interp,ires_dct):
    """
        map hydraulic parameter values to an index array
    """

    hydprop_dict = df_hydprop.to_dict()

    nlay = geol_interp.shape[0]
    nrow = geol_interp.shape[1]
    ncol = geol_interp.shape[2]

    geol_interp_rs = np.reshape(geol_interp,nlay*nrow*ncol)
    geol_array = np.asarray([ires_dct[i] for i in geol_interp_rs])
    par_array = np.zeros((len(hydprop_dict),nlay*nrow*ncol),dtype=float)

    j = 0
    for par, paritems in hydprop_dict.items():
        par_array[j] = [hydprop_dict[par][i] for i in geol_array]
        j=j+1

    par_array = np.reshape(par_array,(len(hydprop_dict),nlay,nrow,ncol))

    return par_array

def write_k_s(sim,par_array,vani=False):

    # Number of layers
    nlay = sim.nlay

    # Model dir
    model_dir = sim.model_ws

    kx_array = par_array[0]
    kz_array = par_array[1]
    ss_array = par_array[2]
    sy_array = par_array[3]

    kx_fnames = []
    kz_fnames = []
    ss_fnames = []
    sy_fnames = []

    for i in range(nlay):
        np.savetxt(os.path.join(model_dir,"kx_"+str(i+1)+".dat"),kx_array[i])
        if vani:
            np.savetxt(os.path.join(model_dir,"vani_"+str(i+1)+".dat"),kz_array[i])
            kz_fnames.append("vani_"+str(i+1)+".dat")
        else:
            np.savetxt(os.path.join(model_dir,"kz_"+str(i+1)+".dat"),kz_array[i])
            kz_fnames.append("kz_"+str(i+1)+".dat")
        np.savetxt(os.path.join(model_dir,"ss_"+str(i+1)+".dat"),ss_array[i])
        np.savetxt(os.path.join(model_dir,"sy_"+str(i+1)+".dat"),sy_array[i])
        kx_fnames.append("kx_"+str(i+1)+".dat")
        ss_fnames.append("ss_"+str(i+1)+".dat")
        sy_fnames.append("sy_"+str(i+1)+".dat")


def eval_geology_east(sim, geol_interp,shape_east,perc=0.3):
    """
        update geology evaluation considering east of Atacama fault
    """
    geol_interp2=geol_interp.copy()
    
    #Grid discretisation
    nlay = sim.nlay
    delc = sim.dis.delc.get_value()
    delr = sim.dis.delr.get_value()

    # Se inicializa la GridIntersect
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)
    # Celdas intersectadas y areas
    shapeRecs = shape_east.shapeRecords()[0]
    shape = shapeRecs.shape
    p = Polygon(shell=shape.points)
    cells = ix.intersect(p)['cellids']
    areas = ix.intersect(p)['areas']

    for k in range(nlay):
        for id,icell in enumerate(cells):
            if areas[id] > perc*delr*delc:
                geol_interp2[k,icell[0],icell[1]]=geol_interp2[k,icell[0],icell[1]]+100

    return geol_interp2


def eval_geology_basamento(sim,shp_basamento,gls_lays,geol_interp,res_dct):

    #Grid discretisation
    delc = sim.dis.delc.get_value()
    delr = sim.dis.delr.get_value()
    lays = np.cumsum(gls_lays)

    # Se inicializa la GridIntersect
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)
    # Celdas intersectadas y areas
    shapeRecs = shp_basamento.shapeRecords()[0]
    shape = shapeRecs.shape
    p = Polygon(shell=shape.points)
    cells = ix.intersect(p)['cellids']

    for k in range(lays[2],lays[3]):
        for id,icell in enumerate(cells):
                if geol_interp[k-1,icell[0],icell[1]] == res_dct['BasamentoHidrogeologico']:
                    pass
                else:
                    geol_interp[k,icell[0],icell[1]] = res_dct['NucleoTransicionProfundo']

    return geol_interp

def df_well_cell(sim,df_series,df_coord):

    model_z = sim.modelgrid.zcellcenters

    # DF con series
    df_series = df_series.drop('Elev_nivel',axis=1)
    df_series.rename(columns={'Borehole':'ID'},inplace=True) 
    df_series.set_index('ID',inplace=True)
    df_series.drop('MP17-074',inplace=True) #Eliminar pozo fuera del dominio

    # DF con coordenadas
    df_coord.drop(df_coord[df_coord['ID'] == 'MP17-074'].index,inplace=True) #Eliminar pozo fuera del dominio

    df_coord['H_pozo'] = df_coord['Cota_Terreno'] - df_coord['Final_Depth'] 

    # Agregar como columna tuplas de celdas
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True) 
    df_coord['CELL_X'] = df_coord.apply(lambda x: ix.intersects(Point(x['ESTE'],x['NORTE']))['cellids'][0][0], axis=1)
    df_coord['CELL_Y'] = df_coord.apply(lambda x: ix.intersects(Point(x['ESTE'],x['NORTE']))['cellids'][0][1], axis=1)
    df_coord['CELL_Z'] = df_coord.apply(lambda x: np.argmin(np.abs(  model_z[:, x['CELL_X'], x['CELL_Y']] - x['H_pozo'] )), axis=1)
    df_coord['CELL'] = df_coord.apply(lambda x: (x['CELL_Z'],x['CELL_X'],x['CELL_Y']), axis=1)

    df_layer = df_coord.drop(['Final_Depth','Cota_Terreno','H_pozo','CELL_X','CELL_Y','CELL'],axis=1)
    df_layer.rename(columns={'CELL_Z':'LAYER'},inplace=True) 

    df_coord = df_coord.drop(['ESTE','NORTE','Final_Depth','Cota_Terreno','H_pozo','CELL_X','CELL_Y','CELL_Z'], axis=1)
    df_coord.set_index('ID',inplace=True)

    # Unir ambos DF 
    df_cell_obs = pd.merge(df_coord, df_series, left_index=True, right_index=True, how='outer')
    df_cell_obs = df_cell_obs.reset_index()

    # Re-indexar con las celdas de los pozos y sus nombres
    df_cell_obs.set_index(['CELL','ID','Fecha'],inplace=True)

    return df_cell_obs,df_layer

def df_well_pest(df_cell_obs):
    # DF a rellenar
    df_new = pd.DataFrame(columns=['CELL', 'ID', 'Fecha', 'Smooth'])
    df_new = df_new.set_index(['CELL', 'ID','Fecha'])

    # Recorrer Dataframe por celda
    for group, fname  in df_cell_obs.groupby(level=[0]):
        num_wel = len(fname.index.get_level_values(1).unique()) # numero de pozos en la celda
        # Pozos con mas de un valor
        if num_wel >= 2:
            list_wel = fname.index.get_level_values(1).unique().to_list() # lista con nombre de los pozos
            list_num = [len(fname.xs((group,list_wel[a]))) for a in range(len(list_wel))] # cantidad de datos por pozo
            fname2 = fname.xs(list_wel[np.argmax(list_num)],level=1,drop_level=False) # DF con el pozo con más datos
            df_new = pd.concat([df_new,fname2])
        else:
            df_new = pd.concat([df_new,fname])  
            
    df_new.reset_index(inplace=True)
    df_new = df_new.drop('CELL',axis=1) # Eliminar celdas
    df_new['Fecha'] = pd.to_datetime(df_new['Fecha'],dayfirst=True ) # Formato de la fecha

    df_new = df_new[['ID','Fecha','Smooth']]

    return df_new


def select_layer(sim,df_coord,df_pozos_obs,geol_interp,res_dct,dic_gid_UH,dic_gid_k):

    df_coord_UH = df_coord.merge(df_pozos_obs,how='left')
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)
    # Filtrar aquellos que no intersectan
    df_coord_UH['intersecta']=True
    for index,row in df_coord_UH.iterrows():
        try:
            ix.intersects(Point(row['ESTE'],row['NORTE']))['cellids'][0][0]
        except:
            df_coord_UH.loc[index,'intersecta']=False
    df_coord_UH=df_coord_UH[df_coord_UH['intersecta']==True]
    df_coord_UH.reset_index(inplace=True)

    # Top y bot de los pozos
    df_coord_UH['top_well'] = df_coord_UH['Cota_Terreno'] - df_coord_UH['Prof_Antepozo']
    df_coord_UH['bot_well'] = df_coord_UH['Cota_Terreno'] - df_coord_UH['Final_Depth']

    # Intersección grilla
    df_coord_UH['CELL_X'] = df_coord_UH.apply(lambda x: ix.intersects(Point(x['ESTE'],x['NORTE']))['cellids'][0][0], axis=1)
    df_coord_UH['CELL_Y'] = df_coord_UH.apply(lambda x: ix.intersects(Point(x['ESTE'],x['NORTE']))['cellids'][0][1], axis=1)

    df_coord_UH['GID1'] = df_coord_UH['Zona'].map(res_dct) #Zona -> GID1 de la evaluacion geologica

    # Interseccion evaluacion geologica Zona y Unidad
    df_coord_UH['eval_geol'] = df_coord_UH.apply(lambda x: geol_interp[:,x['CELL_X'],x['CELL_Y']], axis=1)
    df_coord_UH['in_geol'] = df_coord_UH.apply(lambda x: x['GID1'] in x['eval_geol'], axis=1)
    df_coord_UH['eval_UH'] = df_coord_UH.apply(lambda x: [dic_gid_UH.get(item,item) for item in x['eval_geol']], axis=1)
    df_coord_UH['in_UH'] = df_coord_UH.apply(lambda x: x['UH'] in x['eval_UH'], axis=1)

    # # Intersección profundidad del pozo
    z_bot = sim.modelgrid.botm
    df_coord_UH['eval_depth'] = df_coord_UH.apply(lambda x: z_bot[:,x['CELL_X'],x['CELL_Y']], axis=1) #bot de todas las celdas para un row col seleccionado
    # seleccionar el top de la celda mas cercana dependiendo del signo 
    df_coord_UH['dif_top'] = df_coord_UH.apply(lambda x: x['top_well'] - x['eval_depth'],axis=1)     
    df_coord_UH['min_top'] = df_coord_UH.apply(lambda x: np.argmin(np.abs(x['dif_top'])),axis=1)   
    df_coord_UH['sign_top'] = df_coord_UH.apply(lambda x: np.sign(x['dif_top'][x['min_top']]),axis=1)   
    df_coord_UH['top_layer'] = df_coord_UH.apply(lambda x: x['min_top'] if x['sign_top'] > 0 else x['min_top']+1 , axis=1) 
    # seleccionar el bot de la celda mas cercana dependiendo del signo 
    df_coord_UH['dif_bot'] = df_coord_UH.apply(lambda x: x['bot_well'] - x['eval_depth'],axis=1)     
    df_coord_UH['min_bot'] = df_coord_UH.apply(lambda x: np.argmin(np.abs(x['dif_bot'])),axis=1)   
    df_coord_UH['sign_bot'] = df_coord_UH.apply(lambda x: np.sign(x['dif_bot'][x['min_bot']]),axis=1)   
    df_coord_UH['bot_layer'] = df_coord_UH.apply(lambda x: x['min_bot'] if x['sign_bot'] > 0 else x['min_bot']+1 , axis=1) 

    df_coord_UH['layer_well'] = df_coord_UH.apply(lambda x: np.arange(x['top_layer'],x['bot_layer']+1 ),axis = 1)
    # Intersección layer pozo y UHs 
    df_coord_UH['inter_UH_layer'] = df_coord_UH.apply(lambda x:  x['UH'] in list(map(x['eval_UH'].__getitem__,x['layer_well'])) ,axis = 1)


    # Ranking de la seleccion de pozos
    def select_rank(x):
        if x['in_geol']: # Si se identifica la Zona
            return int(1)
        elif not x['in_geol'] and x['in_UH'] and x['inter_UH_layer']:
            return int(2)
        elif not x['in_geol'] and x['in_UH'] and not x['inter_UH_layer']:
            return int(3)
        else:
            return int(4)

    # Seleccion de layers (se le suma 1 por formato para el PEST)
    def select_layer(x):
        if x['rank'] == 1: 
            layers_inter = np.where(x['GID1'] == x['eval_geol'])[0]
            return layers_inter[len(layers_inter)//2] # Seleccionar la capa media parte entera + 1
        elif x['rank'] == 2:
            lst_int = list(map(x['eval_UH'].__getitem__,x['layer_well'])) #Interseccion layer y UHs
            index_int = [index for (index , item) in enumerate(lst_int) if item == x['UH']] # Index de la interseccion con la UH seleccionada
            layer_int =  x['layer_well'][index_int] #Layer de la intersección
            geol_int = x['eval_geol'][layer_int] #Geología de la intersección 
            k_layer = [dic_gid_k.get(item,item) for item in geol_int] #Valor de k para la geología intersectada
            index_kmax = [index for index, item in enumerate(k_layer) if item == np.max(k_layer)] #Maximos valores de K
            layer_kmax = layer_int[index_kmax] #Layers con los valores maximos de K
            return layer_kmax[len(layer_kmax)//2] # Seleccionar la capa media parte entera de los kmax + 1
        elif x['rank'] == 3:
            return x['eval_UH'].index(x['UH']) #Primer layer donde esta la UH 
        else:
            return x['layer_well'][-1] #Layer mas profundo

    # Rank y seleccion layer
    df_coord_UH['rank'] = df_coord_UH.apply(lambda x: select_rank(x), axis=1)
    df_coord_UH['CELL_Z'] = df_coord_UH.apply(lambda x: select_layer(x), axis=1)

    df_coord_UH['CELL'] = df_coord_UH.apply(lambda x: (x['CELL_Z'],x['CELL_X'],x['CELL_Y']), axis=1)

    df_coord_UH['LAYER'] = df_coord_UH.apply(lambda x: x['CELL_Z'] + 1, axis=1)

    df_cell = df_coord_UH[['ID', 'CELL_X', 'CELL_Y','CELL_Z','CELL']]
    df_coord_layer =  df_coord_UH[['ID','ESTE','NORTE','LAYER']]

    return df_cell,df_coord_layer


def filter_date(year_0, month_0, day_0, year_F, month_F, day_F,df_smooth):

    d1 = datetime.date(year_0,month_0,day_0)
    d2 = datetime.date(year_F,month_F,day_F)

    df_smooth = df_smooth[(df_smooth['Fecha']>=d1.strftime("%Y/%m/%d")) & (df_smooth['Fecha']<=d2.strftime("%Y/%m/%d")) ]

    df_smooth.Fecha = pd.to_datetime(df_smooth.Fecha)
    df_smooth.set_index('Fecha',inplace=True)
    df_smooth.index = df_smooth.index.strftime('%d/%m/%Y %H:%M:%S')
    df_smooth.reset_index(inplace=True)

    df_smooth = df_smooth[['ID','Fecha','Smooth']]

    return df_smooth

def df_shp_rch(sim,shp_path,str):
    shp_obj = shp.Reader(shp_path)
    df = pd.DataFrame({'Shape':shp_obj.shapeRecords()})
    if str == 'nat':
        df['N'] = df.apply(lambda x: x['Shape'].record['id'],axis=1)
        df['Name'] = df.apply(lambda x: x['Shape'].record['Cod_Simple'],axis=1)
    elif str == 'poza':
        df['N'] = df.apply(lambda x: x['Shape'].record['N'],axis=1)
        df['Name'] = df.apply(lambda x: x['Shape'].record['Poza'],axis=1)
        df['A_Name'] = df.apply(lambda x: x['Shape'].record['Area'],axis=1)
    elif str == 'acop':
        df['N'] = df.apply(lambda x: x['Shape'].record['N'],axis=1)
        df['Name'] = df.apply(lambda x: x['Shape'].record['Name'],axis=1)
        df['A_Name'] = df.apply(lambda x: x['Shape'].record['Area'],axis=1)
    df['Polygon'] = df.apply(lambda x: Polygon(shell=x['Shape'].shape.points),axis=1)
    ix = GridIntersect(sim.modelgrid, method="vertex", rtree=True)
    df['Cellids'] = df.apply(lambda x: list(ix.intersect(x['Polygon'])['cellids']),axis=1)
    df['Areas'] = df.apply(lambda x: np.array(ix.intersect(x['Polygon'])['areas']),axis=1)
    if str == 'nat':
        df['A_Name'] = df.apply(lambda x: np.sum(x['Areas']),axis=1)
    df['Porc_Area'] = df.apply(lambda x: x['Areas']/x['A_Name'],axis=1)

    return df

def df_cell_areas(sim,df_nat,df_poza,df_acop):

    nrow = sim.nrow
    ncol = sim.ncol

    df_cell = pd.DataFrame()
    df_cell['Cell'] = [(row,col) for row in range(nrow) for col in range(ncol)]

    def rch_isin(cell,df,str):
        list = df.Cellids.tolist()
        N = df.N.tolist() if str == 'nat' else df.Name.tolist()
        Porc_Area = df.Porc_Area.tolist()
        area = df.Areas.tolist()

        N_lst = []
        A_lst = []
        P_lst = []
        for i,value in enumerate(list):
            if cell in value:
                index = value.index(cell)
                N_lst.append(N[i])
                A_lst.append(area[i][index])
                P_lst.append(Porc_Area[i][index])
        return N_lst,A_lst,P_lst

    df_cell[['Nat','A_Nat','P_Nat']] = df_cell.apply(lambda x: rch_isin(x['Cell'],df_nat,'nat'), axis=1,result_type='expand')
    df_cell[['Poza','A_Poza','P_Poza']] = df_cell.apply(lambda x: rch_isin(x['Cell'],df_poza,'poza'), axis=1,result_type='expand')
    df_cell[['Acop','A_Acop','P_Acop']] = df_cell.apply(lambda x: rch_isin(x['Cell'],df_acop,'acop'), axis=1,result_type='expand')

    # Actualizar valor de area de la recarga natural para que sea todo el area sobrante de la celda
    A = sim.dis.delc.get_value()*sim.dis.delr.get_value()
    df_cell['A_Nat'] = df_cell.apply(lambda x: (A - np.sum(x['A_Poza']) - np.sum(x['A_Acop']))*(x['A_Nat']/np.sum(x['A_Nat'])),axis=1 )

    return df_cell

def conc_SP(SP,inf_pozas,qca_pct):
    
    perlen = SP[0]
    fecha_inicialSP = SP[1]
    fecha_finalSP = SP[2]

    index_df = inf_pozas.index
    df_pond = (inf_pozas.reset_index(drop=True).mul(qca_pct.reset_index(drop=True))).set_index(index_df) #ponderar qca y volumen

    pond_SP = np.zeros((len(perlen),len(df_pond.columns)))
    vol_SP = np.zeros((len(perlen),len(df_pond.columns)))

    for i in range(len(fecha_finalSP)):
        fecha_final = datetime.datetime.strptime(fecha_finalSP[i], "%Y-%m-%d") - datetime.timedelta(days=1) #ARREGLAR
        pond_SP[i,:] = np.sum(df_pond[fecha_inicialSP[i]:fecha_final],axis=0)
        vol_SP[i,:] = np.sum(inf_pozas[fecha_inicialSP[i]:fecha_final],axis=0)

    df_conc_SP = pd.DataFrame(data=np.divide(pond_SP,vol_SP),index=fecha_inicialSP,columns=df_pond.columns) #sum(Ci*Vi)/sum(Vi) por SP
    df_conc_SP.fillna(0.0,inplace=True)
    df_vol_SP = pd.DataFrame(data=vol_SP,index=fecha_inicialSP,columns=df_pond.columns)

    return df_conc_SP,df_vol_SP

def qca_poza(df_cell_rch,df_vol_SP,df_qca_SP,qca):

    def pond_qca_cell(df_vol_SP,df_qca_SP,pozas_cell,pozas_porc):
        if len(pozas_cell)>0:
            vol_porc = df_vol_SP[pozas_cell].reset_index(drop=True).multiply(pozas_porc,axis=1) #volumen ponderado por %area
            cond_pond_cell = (df_qca_SP[pozas_cell].reset_index(drop=True).mul(vol_porc).sum(axis=1)/(vol_porc.sum(axis=1))).fillna(0.0).to_numpy()
            return cond_pond_cell
        else:
            return np.zeros(len(df_vol_SP))

    if qca == 'Li':
        df_cell_rch['Li_Poza'] = df_cell_rch.apply(lambda x: pond_qca_cell(df_vol_SP,df_qca_SP,x['Poza'],x['P_Poza']),axis=1)
    elif qca == 'SO4':
        df_cell_rch['SO4_Poza'] = df_cell_rch.apply(lambda x: pond_qca_cell(df_vol_SP,df_qca_SP,x['Poza'],x['P_Poza']),axis=1)
    elif qca == 'Ca':
        df_cell_rch['Ca_Poza'] = df_cell_rch.apply(lambda x: pond_qca_cell(df_vol_SP,df_qca_SP,x['Poza'],x['P_Poza']),axis=1)
    elif qca == 'K':
        df_cell_rch['K_Poza'] = df_cell_rch.apply(lambda x: pond_qca_cell(df_vol_SP,df_qca_SP,x['Poza'],x['P_Poza']),axis=1)



def qca_array(sim,df_cell_rch):

    nrow = sim.nrow
    ncol = sim.ncol
    nper = sim.nper

    K_pond = np.zeros((nper,nrow,ncol))
    Li_pond = np.zeros((nper,nrow,ncol))
    SO4_pond = np.zeros((nper,nrow,ncol))
    Ca_pond = np.zeros((nper,nrow,ncol))

    df_cell_rch.set_index('Cell',inplace=True)

    for i in range(nrow):
        for j in range(ncol):
            Li_pond[:,i,j] = df_cell_rch.xs((i,j)).Li_Poza
            SO4_pond[:,i,j] = df_cell_rch.xs((i,j)).SO4_Poza
            Ca_pond[:,i,j] = df_cell_rch.xs((i,j)).Ca_Poza
            K_pond[:,i,j] = df_cell_rch.xs((i,j)).K_Poza

    df_cell_rch.reset_index(inplace=True)

    dic_qca = {'Li':Li_pond,'SO4':SO4_pond,'Ca':Ca_pond,'K':K_pond}

    return dic_qca

def chepica_unidad(gls_lays,geol_interp,res_dct):
    u_layer = np.insert(np.cumsum(gls_lays),0,0)
    unidad = ['_UA','_UB','_UC','_UD']

    # Agregar Chepica por Unidad y eliminar el generico
    for i,iuni in enumerate(unidad):
        res_dct.update({'Chepica'+iuni:list(res_dct.values())[-1]+1})
        geol_interp[u_layer[i]:u_layer[i+1],np.where(geol_interp[u_layer[i]:u_layer[i+1],:,:]==6)[1],np.where(geol_interp[u_layer[i]:u_layer[i+1],:,:]==6)[2]] = res_dct['Chepica'+iuni]
    # res_dct.pop('Chepica',None)

    return res_dct,geol_interp