from core.configuration import ConfigurationObject

name = "SAHM"
identifier = "gov.usgs.sahm"
version = '0.0.4'

sahm_path = None
models_path = None
configuration = \
    ConfigurationObject(output_dir= r'I:\VisTrails\WorkingFiles\workspace',
                        layers_path = r'I:\VisTrails\Central_VisTrailsInstall\VisTrails\vistrails\packages\sahm\layers.csv',
                        r_path = r'I:\VisTrails\Central_VisTrailsInstall\Central_R\R-2.12.1\bin',
                        gdal_path = r'I:\VisTrails\Central_VisTrailsInstall\Central_GDAL',
                        verbose = 'True')
