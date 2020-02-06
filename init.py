#  -*- coding: latin-1 -*-
###############################################################################
#  # Contact: talbertc@usgs.gov
#  #
#  # This file is part of the Software for Assisted Habitat Modeling package
#  # for VisTrails.
#  #
#  # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#  # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
#  # THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
#  # PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
#  # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#  # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
#  # OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
#  # WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
#  # OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
#  # ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
#  #
#  # Although this program has been used by the U.S. Geological Survey (USGS),
#  # no warranty, expressed or implied, is made by the USGS or the
#  # U.S. Government as to the accuracy and functioning of the program and
#  # related program material nor shall the fact of distribution constitute
#  # any such warranty, and no responsibility is assumed by the USGS
#  # in connection therewith.
#  #
#  # Any use of trade, firm, or product names is for descriptive purposes only
#  # and does not imply endorsement by the U.S. Government.
###############################################################################

import csv
import os, sys
import shutil
import subprocess
import copy
import time

from vistrails.core.modules.vistrails_module import Module, ModuleError, ModuleSuspended
from vistrails.core.modules.basic_modules import File, Path, Constant
from vistrails.gui.modules.module_configure import StandardModuleConfigurationWidget
from vistrails.core.packagemanager import get_package_manager
import vistrails.core.upgradeworkflow as upgradeworkflow
UpgradeModuleRemap = upgradeworkflow.UpgradeModuleRemap
from vistrails.core import system

from PyQt4 import QtCore, QtGui

from widgets import get_predictor_widget, get_predictor_config

from SelectPredictorsLayers import SelectListDialog
from SelectAndTestFinalModel import SelectAndTestFinalModel

import utils
import GenerateModuleDoc as GenModDoc
import pySAHM.utilities as utilities
import pySAHM.FieldDataAggreagateAndWeight as FDAW
import pySAHM.MDSBuilder as MDSB
import pySAHM.PARC as parc
import pySAHM.RasterFormatConverter as RFC
import pySAHM.SpatialUtilities as SpatialUtilities
from SahmOutputViewer import ModelOutputViewer, ResponseCurveExplorer
from SahmSpatialOutputViewer import ModelMapViewer

from spatial_modules import BaseGeoViewerCell, GeoSpatialViewerCell, RasterLayer, \
                            VectorLayer, PolyLayer, PointLayer

from utils import writetolog
from pySAHM.utilities import TrappedError

identifier = 'gov.usgs.sahm'

doc_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "documentation.xml"))
GenModDoc.load_documentation(doc_file)

def menu_items():
    """ Add a menu item which allows users to specify their session directory
    and select and test the final model
    """
    def change_session_folder():
        global session_dir

        path = str(QtGui.QFileDialog.getExistingDirectory(None,
                                        'Browse to new session folder -', utils.getrootdir()))
        if path == '':
            return None

        session_dir = path
        utils.setrootdir(path)
        utils.createLogger(session_dir, True)

        package_manager = get_package_manager()
        package = package_manager.get_package(identifier)
        configuration.set_deep_value('cur_session_folder', path)
        package.persist_configuration()

        writetolog("*" * 79 + "\n" + "*" * 79)
        writetolog(" output directory:   " + session_dir)
        writetolog("*" * 79 + "\n" + "*" * 79)

    def select_test_final_model():
        global session_dir

        STFM = SelectAndTestFinalModel(session_dir, utils.get_r_path())
        retVal = STFM.exec_()

    def selectProcessingMode():
        selectDialog = QtGui.QDialog()

        global groupBox
        groupBox = QtGui.QGroupBox("Processing mode:")
        vbox = QtGui.QVBoxLayout()

        for mode in [("multiple models simultaneously (1 core each)", True),
                     ("single models sequentially (n - 1 cores each)", True)]:
            radio = QtGui.QRadioButton(mode[0])
            radio.setChecked(mode[0] == configuration.cur_processing_mode)
            radio.setEnabled(mode[1])
            QtCore.QObject.connect(radio, QtCore.SIGNAL("toggled(bool)"), selectProcessingMode_changed)
            vbox.addWidget(radio)

        groupBox.setLayout(vbox)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(groupBox)
        selectDialog.setLayout(layout)

        selectDialog.exec_()

    def selectProcessingMode_changed(e):
        if e:
            global groupBox
            qvbl = groupBox.layout()
            for i in range(0, qvbl.count()):
                widget = qvbl.itemAt(i).widget()
                if (widget != 0) and (type(widget) is QtGui.QRadioButton):
                    if widget.isChecked():
                        package_manager = get_package_manager()
                        package = package_manager.get_package(identifier)
                        configuration.set_deep_value('cur_processing_mode', str(widget.text()))
                        package.persist_configuration()

                        utilities.start_new_pool(utilities.get_process_count(widget.text()))


    lst = []
    lst.append(("Change session folder", change_session_folder))
    lst.append(("Change processing mode", selectProcessingMode))
    lst.append(("Select and test the Final Model", select_test_final_model))
    return(lst)

class SAHMDocumentedModule(object):

    @classmethod
    def provide_input_port_documentation(cls, port_name):
        return GenModDoc.construct_port_doc(cls, port_name, 'in')
    @classmethod
    def provide_output_port_documentation(cls, port_name):
        return GenModDoc.construct_port_doc(cls, port_name, 'out')

class SAHMPathModule(Module, SAHMDocumentedModule):
    '''The base class that all Path modules in SAHM inherit from
    '''
    _input_ports = [('file', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}), ]

    def compute(self):
        self.set_output('value', utils.get_relative_path(self.get_input("file"), module=self))

class FieldData(SAHMPathModule):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('FieldData')
    _output_ports = [('value', '(gov.usgs.sahm:FieldData:DataInput)'), ]


class TemplateLayer(Module, SAHMDocumentedModule):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('TemplateLayer')
    _input_ports = [('file', '(edu.utah.sci.vistrails.basic:File)', {'optional':True})]
    _output_ports = [('value', '(gov.usgs.sahm:TemplateLayer:DataInput)')]
    
    def compute(self):
        items = self.force_get_input_list('file')
        for i in items:
            if len(str(i)) > 14:
                template = utils.get_relative_path(i, module=self)
                writetolog(str(template))
                self.set_output('value', template)

class Predictor(Module, SAHMDocumentedModule):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('Predictor')

    _input_ports = [('categorical', '(edu.utah.sci.vistrails.basic:Boolean)', {'optional': True}),
                    ('ResampleMethod', '(edu.utah.sci.vistrails.basic:String)',
                     {'entry_types': "['enum']",
                      'values': "[['NearestNeighbor', 'Bilinear', 'Cubic', 'CubicSpline', 'Lanczos']]",
                      'defaults':'["Bilinear"]', 'optional': True}),
                    ('AggregationMethod', '(edu.utah.sci.vistrails.basic:String)',
                     {'entry_types': "['enum']",
                      'values': "[['Mean', 'Max', 'Min', 'STD', 'Majority', 'None']]", 'optional': True,
                      'defaults':'["Mean"]', 'optional': True}),
                    ('file', '(edu.utah.sci.vistrails.basic:Path)', {'optional': True})]
    _output_ports = [('value', '(gov.usgs.sahm:Predictor:DataInput)'), ]

    def compute(self):
        if (self.has_input("ResampleMethod")):
            resample_method = self.get_input("ResampleMethod")
            if resample_method.lower() not in ['nearestneighbor', 'bilinear', 'cubic', 'cubicspline', 'lanczos']:
                raise ModuleError(self,
                                  "Resample Method not one of 'nearestneighbor', 'bilinear', 'cubic', 'cubicspline', or 'lanczos'")
        else:
            resample_method = 'Bilinear'

        if (self.has_input("AggregationMethod")):
            aggregation_method = self.get_input("AggregationMethod")
            if self.get_input("AggregationMethod").lower() not in ['mean', 'max', 'min', 'std', 'majority', 'none']:
                raise ModuleError(self, "No Aggregation Method specified")
        else:
            aggregation_method = "Mean"

        if (self.has_input("categorical")):
            if self.get_input("categorical") == True:
                categorical = '1'
            else:
                categorical = '0'
        else:
            categorical = '0'

        if (self.has_input("file")):
            out_fname = utils.get_relative_path(self.get_input("file"), self)
            inFile = utils.get_raster_name(out_fname)
        else:
            raise ModuleError(self, "No input file specified")
        self.set_output('value', (inFile, categorical, resample_method, aggregation_method))

class PredictorList(Constant):
    '''
    This module is a required class for other modules and scripts within the
    SAHM package. It is not intended for direct use or incorporation into
    the VisTrails workflow by the user.
    '''
    _input_ports = [('value', '(gov.usgs.sahm:PredictorList:Other)'),
                                 ('addPredictor', '(gov.usgs.sahm:Predictor:DataInput)')]
    _output_ports = [('value', '(gov.usgs.sahm:PredictorList:Other)')]

    @staticmethod
    def translate_to_string(v):
        return str(v)

    @staticmethod
    def translate_to_python(v):
        v_list = eval(v)
        return v_list

    @staticmethod
    def validate(x):
        return type(x) == list

    def compute(self):
        p_list = self.force_get_input_list("addPredictor")
        v = self.force_get_input("value", [])

        b = self.validate(v)
        if not b:
            raise ModuleError(self, "Internal Error: Constant failed validation")
        if len(v) > 0 and type(v[0]) == tuple:
            f_list = [utils.get_relative_path(v_elt[0], module=self) for v_elt in v]
        else:
            f_list = v
        p_list += f_list
        #  self.set_output("value", p_list)
        self.set_output("value", v)

class PredictorListFile(SAHMDocumentedModule, Module):
    '''
    copies the input predictor list csv to our working directory
    and appends any additionally added predictors
    '''
    __doc__ = GenModDoc.construct_module_doc('PredictorListFile')

    _input_ports = [('csvFileList', '(edu.utah.sci.vistrails.basic:File)', {'optional': True})]
    _output_ports = [('RastersWithPARCInfoCSV', '(gov.usgs.sahm:RastersWithPARCInfoCSV:Other)')]

    def compute(self):
        if not self.has_input("csvFileList"):
            raise ModuleError(self, "No CSV file provided")

        output_file = utils.get_relative_path(self.get_input("csvFileList"), module=self)
        self.set_output('RastersWithPARCInfoCSV', output_file)





#  ##Internal class definitions, used for port connection enforcement

class MergedDataSet(File):
    '''
    This module is a required class for other modules and scripts within the
    SAHM package. It is not intended for direct use or incorporation into
    the VisTrails workflow by the user.
    '''
    _input_ports = [('mdsFile', '(edu.utah.sci.vistrails.basic:File)'), ]
    _output_ports = [('value', '(gov.usgs.sahm:MergedDataSet:Other)'), ]

    pass

class RastersWithPARCInfoCSV(File):
    '''
    This module is a required class for other modules and scripts within the
    SAHM package. It is not intended for direct use or incorporation into
    the VisTrails workflow by the user.
    '''
    _input_ports = [('mdsFile', '(edu.utah.sci.vistrails.basic:File)'), ]
    _output_ports = [('value', '(gov.usgs.sahm:MergedDataSet:Other)'), ]

    pass


#  ##End of Internal class definitions, used for port connection enforcement


class Model(SAHMDocumentedModule, Module):
    '''
    This module is a required class for other modules and scripts within the
    SAHM package. It is not intended for direct use or incorporation into
    the VisTrails workflow by the user.
    '''
    _input_ports = [('ThresholdOptimizationMethod', '(edu.utah.sci.vistrails.basic:String)',
                     {'entry_types': "['enum']",
                      'values': '[["Threshold=0.5", "Sensitivity=Specificity", "Maximizes (sensitivity+specificity)/2", "Maximizes Cohen\'s Kappa","Maximizes PCC (percent correctly classified)","Predicted prevalence=observed prevalence","Threshold=observed prevalence","Mean predicted probability","Minimizes distance between ROC plot and (0,1)",]]', 'optional': True,
                      'defaults':'["Sensitivity=Specificity"]', 'optional':True}),
                    ('mdsFile', '(gov.usgs.sahm:MergedDataSet:Other)'),
                    ('makeBinMap', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                    ('makeProbabilityMap', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                    ('makeMESMap', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["False"]', 'optional':True}),
                    ('outputFolderName', '(edu.utah.sci.vistrails.basic:String)', {'optional':True}),
                    ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}), ]

    _output_ports = [('modelWorkspace', '(edu.utah.sci.vistrails.basic:Directory)'),
                     ('BinaryMap', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                     ('ProbabilityMap', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                     ('ResidualsMap', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                     ('MessMap', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                     ('MoDMap', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                     ('modelEvalPlot', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                     ('Text_Output', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                     ('ModelVariableImportance', '(edu.utah.sci.vistrails.basic:File)', {'optional': True}),]

    port_map = {'mdsFile':('c', None, True),  #  These ports are for all Models
                         'makeProbabilityMap':('mpt', utils.R_boolean, True),
                         'makeBinMap':('mbt', utils.R_boolean, True),
                         'makeMESMap':('mes', utils.R_boolean, True),
                         'ThresholdOptimizationMethod':('om', None, False),
                    }

    def __init__(self):
        self.suspended_completed = False
        self.pywrapper = "runRModel.py"
        self.port_map = copy.deepcopy(Model.port_map)
        self.output_dname = None
        Module.__init__(self)

    def compute(self):
        out_folder = self.force_get_input("outputFolderName", "")

        self.args_dict = utils.map_ports(self, self.port_map)

        mdsFile = utils.get_relative_path(self.args_dict['c'], self)

        if self.has_input('run_name_info'):
            runinfo = self.force_get_input('run_name_info')
            if not type(runinfo) == dict:
                runinfo = runinfo.contents
            subfolder = runinfo.get('subfolder_name', "")
            runname = runinfo.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(mdsFile)

        prefix = self.abbrev
        if prefix == "ApplyModel":
            prefix = "Apply"
            if not os.path.isdir(self.args_dict['ws']):
                orig_ws = os.path.split(self.args_dict['ws'])[0]
            else:
                orig_ws = self.args_dict['ws']
            prefix += os.path.split(orig_ws)[1].replace('_', '').upper()

        name_items = filter(None, [prefix, runname, out_folder])
        prefix = "_".join(name_items)

        #  convert threshold optimization string to the expected integer
        thresholds = {"Threshold=0.5":1,
                      "Sensitivity=Specificity":2,
                      "Maximizes (sensitivity+specificity)/2":3,
                      "Maximizes Cohen's Kappa":4,
                      "Maximizes PCC (percent correctly classified)":5,
                      "Predicted prevalence=observed prevalence":6,
                      "Threshold=observed prevalence":7,
                      "Mean predicted probability":8,
                      "Minimizes distance between ROC plot and (0,1)":9}
        self.args_dict["om"] = thresholds.get(self.args_dict.get("om", "Sensitivity=Specificity"))

        if not utils.checkModelCovariatenames(mdsFile):
            msg = "These R models do not work with covariate names begining with non-letter characters or \n"
            msg += "covariate names that contain non-alphanumeric characters other than '.' or '_'.\n"
            msg += "Please check your covariate names and rename any that do not meet these constraints.\n"
            msg += "Covaraiate names are found in the first line of the mds file: \n"
            msg += "\t\t" + mdsFile
            writetolog(msg, False, True)
            raise ModuleError(self, msg)

        if self.abbrev == "Maxent":
            self.args_dict['maxent_path'] = configuration.maxent_path
            self.args_dict['java_path'] = utils.find_java_exe(configuration.java_path)
            self.args_dict['maxent_args'] = self.maxent_args

        self.args_dict['rc'] = utils.MDSresponseCol(mdsFile)
        self.args_dict['cur_processing_mode'] = configuration.cur_processing_mode


        self.output_dname, signature, already_run = utils.make_next_file_complex(self, prefix, key_inputs=[mdsFile],
                                                                                file_or_dir='dir', subfolder=subfolder)
        copy_mds_fname = os.path.join(self.output_dname, os.path.split(mdsFile)[1])
        if not os.path.exists(copy_mds_fname):
            shutil.copyfile(mdsFile, copy_mds_fname)
        expanded_output_dname = os.path.join(self.output_dname, "ExpandedOutput")
        if not os.path.exists(expanded_output_dname):
            os.makedirs(expanded_output_dname)

        self.args_dict["c"] = copy_mds_fname

#            self.output_dname = utils.find_model_dir(prefix, self.args_dict)

        if self.abbrev == 'brt' or \
            self.abbrev == 'rf':
            if not "seed" in self.args_dict.keys():
                self.args_dict['seed'] = utils.get_seed()
            writetolog("    seed used for " + self.abbrev + " = " + str(self.args_dict['seed']))

        self.args_dict['o'] = self.output_dname

        #  This give previously launched models time to finish writing their
        #  logs so we don't get a lock
        time.sleep(0.5)

        utils.write_hash_entry_pickle(signature, self.output_dname)

        try:
            utils.run_model_script(self.name, self.args_dict, self, self.pywrapper)
        except ModuleSuspended:
            raise
        except:
            utils.delete_hash_entry_pickle(signature)
            raise

        self.set_model_results()

    def set_model_results(self,):
        #  set our output ports
        #  if an output is expected and we're running in syncronously then throw
        #  an error
        if not self.args_dict.has_key('mes'):
            self.args_dict['mes'] = 'FALSE'
        self.outputRequired = configuration.cur_processing_mode == "single models sequentially (n - 1 cores each)"
        self.setModelResult("_prob_map.tif", 'ProbabilityMap', self.abbrev)
        self.setModelResult("_bin_map.tif", 'BinaryMap', self.abbrev)
        self.setModelResult("_resid_map.tif", 'ResidualsMap', self.abbrev)
        self.setModelResult("_mess_map.tif", 'MessMap', self.abbrev)
        self.setModelResult("_MoD_map.tif", 'MoDMap', self.abbrev)
        self.setModelResult("_output.txt", 'Text_Output', self.abbrev)
        self.setModelResult("_modelEvalPlot.png", 'modelEvalPlot', self.abbrev)
        self.setModelResult("_variable.importance.png", 'ModelVariableImportance', self.abbrev)
        writetolog("Finished " + self.abbrev + " builder\n", True, True)

        self.set_output("modelWorkspace", self.output_dname)

    def setModelResult(self, filename, portname, abbrev):
        '''sets a single output port value
        '''
        out_fname = os.path.join(self.output_dname, abbrev + filename)
        self.set_output(portname, out_fname)

class GLM(Model):
    __doc__ = GenModDoc.construct_module_doc('GLM')

    _input_ports = list(Model._input_ports)
    _input_ports.extend([('SelectBestPredSubset', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                         ('SimplificationMethod', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["AIC"]', 'optional':True}),
                         ('SquaredTerms', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                         ])
    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'FIT_GLM_pluggable.r'
        self.abbrev = 'glm'
        self.port_map.update({'SimplificationMethod':('sm', None, False),  #  This is a GLM specific port
                         'SquaredTerms':('sqt', utils.R_boolean, False),  #  This is a GLM specific port
                         'SelectBestPredSubset':('pst', utils.R_boolean, False),  #  This is a GLM specific port
                         })

class RandomForest(Model):
    __doc__ = GenModDoc.construct_module_doc('RandomForest')

    _input_ports = list(Model._input_ports)
    _input_ports.extend([('Seed', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["{}"]'.format(utils.get_seed()), 'optional':True}),
                         ('mTry', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["1"]', 'optional':True}),
                         ('nTrees', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                         ('nodesize', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                         ('replace', '(edu.utah.sci.vistrails.basic:Boolean)', {'optional':True}),
                         ('maxnodes', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                         ('importance', '(edu.utah.sci.vistrails.basic:Boolean)', {'optional':True}),
                         ('localImp', '(edu.utah.sci.vistrails.basic:Boolean)', {'optional':True}),
                         ('proximity', '(edu.utah.sci.vistrails.basic:Boolean)', {'optional':True}),
                         ('oobProx', '(edu.utah.sci.vistrails.basic:Boolean)', {'optional':True}),
                         ('normVotes', '(edu.utah.sci.vistrails.basic:Boolean)', {'optional':True}),
                         ])
    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'FIT_RF_pluggable.r'
        self.abbrev = 'rf'
        self.port_map.update({'Seed':('seed', utils.get_seed, True),  #  This is a BRT specific port
                         'mTry': ('mtry', None, False),  #  This is a Random Forest specific port
                         'nTrees': ('ntree', None, False), # MMF 02/22/2016 - nTrees was not getting passed
                         'nodesize': ('nodeS', None, False),  #  This is a Random Forest specific port
                         'replace': ('sampR', utils.R_boolean, False),  #  This is a Random Forest specific port
                         'maxnodes': ('maxN', None, False),  #  This is a Random Forest specific port
                         'importance': ('impt', utils.R_boolean, False),  #  This is a Random Forest specific port
                         'localImp': ('locImp', utils.R_boolean, False),  #  This is a Random Forest specific port
                         'proximity': ('prox', utils.R_boolean, False),  #  This is a Random Forest specific port
                         'oobPorx': ('oopp', utils.R_boolean, False),  #  This is a Random Forest specific port
                         'normVotes': ('nVot', utils.R_boolean, False),  #  This is a Random Forest specific port
                         'doTrace': ('Trce', utils.R_boolean, False),  #  This is a Random Forest specific port
                         'keepForest': ('kf', utils.R_boolean, False),  #  This is a Random Forest specific port
                         })

class MARS(Model):
    __doc__ = GenModDoc.construct_module_doc('MARS')

    _input_ports = list(Model._input_ports)
    _input_ports.extend([('MarsDegree', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["1"]', 'optional':True}),
                         ('MarsPenalty', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["2.0"]', 'optional':True}),
                          ])
    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'FIT_MARS_pluggable.r'
        self.abbrev = 'mars'
        self.port_map.update({'MarsDegree':('deg', None, False),  #  This is a MARS specific port
                         'MarsPenalty':('pen', None, False),  #  This is a MARS specific port
                         })

class ApplyModel(Model):
    __doc__ = GenModDoc.construct_module_doc('ApplyModel')
    _input_ports = list(Model._input_ports)
    _input_ports.insert(0, ('modelWorkspace', '(edu.utah.sci.vistrails.basic:Directory)'))
#    _input_ports.insert(1, ('evaluateHoldout', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["False"]', 'optional':False}))
#    _input_ports.extend([('modelWorkspace', '(edu.utah.sci.vistrails.basic:Directory)')])

    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'EvaluateNewData.r'
        self.abbrev = 'ApplyModel'
        self.port_map.update({'modelWorkspace':('ws',
                lambda x: os.path.join(utils.dir_path_value(x), "modelWorkspace"), True), })

    def compute(self):
        #  if the suplied mds has rows, observations then
        #  pass r code the flag to produce metrics
        mdsfname = utils.get_relative_path(self.force_get_input('mdsFile'), self)
        workspace = utils.get_relative_path(self.force_get_input('modelWorkspace'), self)

        mdsfile = open(mdsfname, "r")
        lines = mdsfile.readlines()

        if len(lines) > 3:
            #  we have rows R will need to recreate metrics.
            self.args = 'pmt=TRUE '
        else:
            self.args = 'pmt=FALSE '

        #  make sure all the covariates in the original model are in the new csv
        #  if not raise an exception alerting the user
        orig_mds = utils.get_mdsfname(workspace)
        orig_mdsfile = open(orig_mds, "r")
        orig_lines = orig_mdsfile.readlines()


        skip_list = ['Split', 'EvalSplit', 'Weights']

        orig_covariates = [item.strip() for item in orig_lines[0].split(",")[3:]
                                            if item.strip() not in skip_list]
        orig_used = [item.strip() for item in orig_lines[1].split(",")[3:]]
        missing_covariates = []
        new_covariates = [item.strip() for item in lines[0].split(",")[3:]
                                            if item.strip() not in skip_list]

        for orig_covariate in orig_covariates:
            i = orig_covariates.index(orig_covariate)
            if orig_used[i] == "1" and \
                    new_covariates.count(orig_covariate) == 0:
                missing_covariates.append(orig_covariate)

        if len(missing_covariates) > 0:
            msg = 'One or more of the covariates used in the original model are not specified in the apply model mds file\n'
            msg += 'Specfically the following covariates were not found:'
            msg += '\n\t'.join(missing_covariates)

            raise RuntimeError(msg)

        Model.compute(self)

class xgBoost(Model):
    __doc__ = GenModDoc.construct_module_doc('xgBoost')

    _input_ports = list(Model._input_ports)
    _input_ports.extend([('nrounds', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["10"]', 'optional':True}),
                         ('eta', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["1.0"]', 'optional':True}),
                         ('gamma', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["0"]', 'optional':True}),
                         ('max_depth', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["6"]', 'optional':True}),                         
                         ('min_child_weight', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["0"]', 'optional':True}),
                         ('subsample', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["0.632"]', 'optional':True}),
                          ])

    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'FIT_XGB_pluggable.r'
        self.abbrev = 'xgb'
        self.port_map.update({'nrounds':('nrounds', None, False),  #  These are all xgBoost specific ports
                         'eta':('eta', None, False),               
                         'gamma':('gamma', None, False),           
                         'max_depth':('md', None, False),          
                         'min_child_weight':('mcw', None, False),  
                         'subsample':('samp', None, False),        
                         })

class BoostedRegressionTree(Model):
    __doc__ = GenModDoc.construct_module_doc('BoostedRegressionTree')

    _input_ports = list(Model._input_ports)
    _input_ports.extend([('Seed', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["{}"]'.format(utils.get_seed()), 'optional':True}),
                              ('TreeComplexity', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                              ('BagFraction', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["0.75"]', 'optional':True}),
                              ('NumberOfFolds', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["3"]', 'optional':True}),
                              ('Alpha', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["1"]', 'optional':True}),
                              ('PrevalenceStratify', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                              ('ToleranceMethod', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["auto"]', 'optional':True}),
                              ('Tolerance', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["0.001"]', 'optional':True}),
                              ('LearningRate', '(edu.utah.sci.vistrails.basic:Float)', {'optional':True}),
                              ('SelectBestPredSubset', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                              ('NumberOfTrees', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                              ])
    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'FIT_BRT_pluggable.r'
        self.abbrev = 'brt'
        self.port_map.update({'Seed':('seed', None, False),  #  This is a BRT specific port
                         'TreeComplexity':('tc', None, False),  #  This is a BRT specific port
                         'BagFraction':('bf', None, False),  #  This is a BRT specific port
                         'NumberOfFolds':('nf', None, False),  #  This is a BRT specific port
                         'Alpha':('alp', None, False),  #  This is a BRT specific port
                         'PrevalenceStratify':('ps', None, False),  #  This is a BRT specific port
                         'ToleranceMethod':('tolm', None, False),  #  This is a BRT specific port
                         'Tolerance':('tol', None, False),  #  This is a BRT specific port
                         'LearningRate':('lr', None, False),  #  This is a BRT specific port
                         'NumberOfTrees':('ntr', None, False),  #  This is a BRT specific port
                         'SelectBestPredSubset':('pst', utils.R_boolean, False),  #  This is a BRT specific port
                         })

class MAXENT(Model):
    '''
    '''
    _input_ports = list(Model._input_ports)
    _input_ports.extend([('UseRMetrics', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                         ])
    _output_ports = list(Model._output_ports)
    _output_ports.extend([("lambdas", "(edu.utah.sci.vistrails.basic:File)", {'optional':True}),
                     ("report", "(edu.utah.sci.vistrails.basic:File)", {'optional':True}),
                     ("roc", "(edu.utah.sci.vistrails.basic:File)", {'optional':True}),])

    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'WrapMaxent.r'
        self.pywrapper = "runMaxent.py"
        self.abbrev = 'Maxent'
        self.port_map.update({'species_name':('species_name', None, True),  #  This is a Maxent specific port
                              })

    def compute(self):

        self.maxent_args = {}
        for port in self._input_ports:

            if not port in list(Model._input_ports) and \
                port[0] <> 'projectionlayers' and \
                port[0] <> 'UseRMetrics' and \
                port[0] <> 'species_name':
                if self.has_input(port[0]):
                    port_val = self.get_input(port[0])
                    if port[1] == "(edu.utah.sci.vistrails.basic:Boolean)":
                        port_val = str(port_val).lower()
                    elif (port[1] == "(edu.utah.sci.vistrails.basic:Path)" or \
                        port[1] == "(edu.utah.sci.vistrails.basic:File)" or \
                        port[1] == "(edu.utah.sci.vistrails.basic:Directory)"):
                        port_val = port_val.name
                    self.maxent_args[port[0]] = port_val
                else:
                    kwargs = port[2]
                    try:
                        if port[1] == "(edu.utah.sci.vistrails.basic:Boolean)":
                            default = kwargs['defaults'][1:-1].lower()
                        elif port[1] == "(edu.utah.sci.vistrails.basic:String)":
                            default = kwargs['defaults'][1:-1]
                        else:
                            default = kwargs['defaults'][1:-1]
                        #  args[port[0]] = default
                        self.maxent_args[port[0]] = default[1:-1]
                    except KeyError:
                        pass
        if self.has_input('projectionlayers'):
            value = self.force_get_input_list('projectionlayers')
            projlayers = ','.join([path.name for path in value])
            self.maxent_args['projectionlayers'] = projlayers

        Model.compute(self)

#       set some Maxent specific outputs
        self.args_dict['species_name'] = self.args_dict['species_name'].replace(' ', '_')
        lambdasfile = self.args_dict["species_name"] + ".lambdas"
        self.setModelResult(lambdasfile, "lambdas", "")

        rocfile = "plots" + os.sep + self.args_dict["species_name"] + "_roc.png"
        self.setModelResult(rocfile, "roc", "")

        htmlfile = self.args_dict["species_name"] + ".html"
        self.setModelResult(htmlfile, "report", "")

        writetolog("Finished Maxent", True)

class UserDefinedCurve(Model):
    __doc__ = GenModDoc.construct_module_doc('UserDefinedCurve')

    _input_ports = list(Model._input_ports)
    _input_ports.extend([("curves_json", "(edu.utah.sci.vistrails.basic:File)", {'optional':True}),])
    _output_ports = list(Model._output_ports)
    _output_ports.extend([("curves_json", "(edu.utah.sci.vistrails.basic:File)", {'optional':True}),])

    def __init__(self):
        global models_path
        Model.__init__(self)
        self.name = 'FIT_UDC.r'
        self.pywrapper = "runRModel.py"
        self.abbrev = 'udc'
        
        self.port_map.update({'curves_json':('curves_json', None, False),  #  This is a Maxent specific port
                              })

    def compute(self):

        Model.compute(self)

        self.setModelResult("udc.json", "curves_json", "")

        writetolog("Finished UserDefinedCurves", True)

class EnsembleBuilder(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('EnsembleBuilder')

    _input_ports = [('ModelWorkspaces', '(edu.utah.sci.vistrails.basic:Directory)'),
                    ('ThresholdMetric', '(edu.utah.sci.vistrails.basic:String)',
                     {'entry_types': "['enum']",
                      'values': '[["None", "AUC", "Percent Correctly Classified", "Sensitivity", "Specificity", "Kappa", "True Skill Statistic"]]',
                      'defaults':'["None"]', 'optional':True}),
                    ('ThresholdValue', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["0.75"]', 'optional':True}),
                    ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':True}), ]
    _output_ports = [("AverageProbability", "(edu.utah.sci.vistrails.basic:File)"),
                     ("BinaryCount", "(edu.utah.sci.vistrails.basic:File)"),]

    def compute(self):
        port_map = {'ThresholdMetric': ('ThresholdMetric', None, True),
                    'ThresholdValue': ('ThresholdValue', None, False),
            'run_name_info': ('run_name_info', None, False), }

        params = utils.map_ports(self, port_map)

        model_workspaces = self.get_input_list("ModelWorkspaces")
        if len(model_workspaces) < 2:
            raise RuntimeError('2 or more ModelWorkspaces must be supplied!')

        #  TODO add in check to make sure all models finished successfully
        #  if still running raise module suspended
        workspaces = []
        for model_workspace in model_workspaces:


            rel_workspace = utils.get_relative_path(model_workspace, self)
            if params['ThresholdMetric'] != 'None':
                model_results = utils.get_model_results(rel_workspace)
                param_key = params['ThresholdMetric'].replace(' ', '').lower()
                if float(model_results[param_key]) >= float(params['ThresholdValue']):
                    workspaces.append(os.path.normpath(rel_workspace))
                else:
                    msg = "Model below threshold, Removed from ensemble! :\n"
                    msg += os.path.normpath(rel_workspace)
                    msg += "\n model {} value of {} below threshold of {}".format(params['ThresholdMetric'], model_results[param_key], params['ThresholdValue'])

                    writetolog(msg, True)
            else:
                workspaces.append(os.path.normpath(rel_workspace))

        run_name_info = params.get('run_name_info')
        if run_name_info:
            if not type(run_name_info) == dict:
                run_name_info = run_name_info.contents
            subfolder = run_name_info.get('subfolder_name', "")
            runname = run_name_info.get('runname', "")
        else:
            #  If all subfolders are the same we'll but the model output in the same subfolder
            subfolders = []
            for model_workspace in model_workspaces:
                subfolder, runname = utils.get_previous_run_info(utils.get_relative_path(model_workspace))
                subfolders.append(subfolder)
            if all(x == subfolders[0] for x in subfolders):
                subfolder = subfolders[0]
            else:
                subfolder = ''
            runname = ''

        prefix = "ensemble_prob"
        suffix = ".tif"

        prob_tifs = [os.path.join(ws, utils.find_file(ws, '_prob_map.tif')) for ws in workspaces]
        bin_tifs = [os.path.join(ws, utils.find_file(ws, '_bin_map.tif')) for ws in workspaces]

        output_fname, signature, already_run = utils.make_next_file_complex(self,
                                        prefix=prefix, suffix=suffix,
                                        key_inputs=prob_tifs,
                                        subfolder=subfolder, runname=runname)
        output_fname_bin = output_fname.replace("ensemble_prob", "ensemble_count")

        if already_run:
            writetolog("No change in inputs or parameters using previous run of EnsembleBuilder", True)
        else:

            SpatialUtilities.average_geotifs(prob_tifs, output_fname, None, False, SpatialUtilities.average_nparrays)
            SpatialUtilities.average_geotifs(bin_tifs, output_fname_bin, None, False, SpatialUtilities.sum_nparrays)

        if os.path.exists(utils.get_relative_path(output_fname, self)):
            writetolog("Finished Ensemble generation ", True)
        else:
            msg = "Problem encountered building ensemble maps.  Expected output file not found."
            writetolog(msg, False)
            raise ModuleError(self, msg)

        utils.write_hash_entry_pickle(signature, output_fname)
        self.set_output("AverageProbability", utils.get_relative_path(output_fname, self))
        self.set_output("BinaryCount", utils.get_relative_path(output_fname_bin, self))

class BackgroundSurfaceGenerator(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('BackgroundSurfaceGenerator')

    _input_ports = [('templateLayer', '(gov.usgs.sahm:TemplateLayer:DataInput)'),
                    ('fieldData', '(gov.usgs.sahm:FieldData:DataInput)'),
                        ('method', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["KDE"]', 'optional':True}),
                        ('bandwidthOptimizationMethod', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["adhoc"]', 'optional':True}),
                        ('isopleth', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["95"]', 'optional':True}),
                        ('continuous', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["False"]', 'optional':True}),
                        ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}), ]
    _output_ports = [("KDE", "(edu.utah.sci.vistrails.basic:File)")]
    # _output_ports = [('KDE', '(gov.usgs.sahm:BackgroundSurfaceGenerator:Tools)')]

    def compute(self):
        port_map = {'templateLayer': ('templatefName', None, True),
                    'fieldData': ('fieldData', None, False),
            'method': ('method', None, True),
            'bandwidthOptimizationMethod': ('bandOptMeth', None, True),
            'isopleth': ('isopleth', None, True),
            'continuous': ('continuous', utils.R_boolean, True),
            'run_name_info': ('run_name_info', None, False), }

        kde_params = utils.map_ports(self, port_map)

        run_name_info = kde_params.get('run_name_info')
        if run_name_info:
            if not type(runinfo) == dict:
                runinfo = runinfo.contents
            subfolder = run_name_info.get('subfolder_name', "")
            runname = run_name_info.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(kde_params['fieldData'])

        global models_path
        prefix = os.path.splitext(os.path.split(kde_params["fieldData"])[1])[0]
        suffix = kde_params["method"]
        if kde_params["method"] == "KDE":
            suffix += kde_params["bandOptMeth"]
            if kde_params["continuous"] == "TRUE":
                suffix += "_continuous"
            else:
                suffix += "_iso" + str(kde_params["isopleth"])
        suffix += ".tif"

        output_fname, signature, already_run = utils.make_next_file_complex(self,
                                        prefix=prefix, suffix=suffix,
                                        key_inputs=[kde_params['fieldData'], utils.get_raster_files(kde_params['templatefName'])],
                                        subfolder=subfolder, runname=runname)

        if already_run:
            writetolog("No change in inputs or parameters using previous run of BackgroundSurfaceGenerator", True)
        else:
            args = {"tmplt":kde_params["templatefName"],
                    "i":kde_params["fieldData"],
                    "o":output_fname,
                    "mth":kde_params["method"],
                    "bwopt":kde_params["bandOptMeth"],
                    "ispt":str(kde_params["isopleth"]),
                    "continuous":kde_params["continuous"]}

            utils.run_R_script("PseudoAbs.r", args, self, new_r_path=configuration.r_path)

        if os.path.exists(output_fname):
            output_file = utils.get_relative_path(output_fname, module=self)
            writetolog("Finished KDE generation ", True)
        else:
            msg = "Problem encountered generating KDE.  Expected output file not found."
            writetolog(msg, False)
            raise ModuleError(self, msg)

        utils.write_hash_entry_pickle(signature, output_fname)
        self.set_output("KDE", output_file)

class OutputNameInfo(Constant):
    contents = {}

    @staticmethod
    def translate_to_python(x):
        try:
            runinfo = OutputNameInfo()
            runinfo.contents = {'runname':'',
                                'subfolder_name':str(x),
                                'delete_previous':False}
            return runinfo
        except:
            return None

class OutputName(SAHMDocumentedModule, Module):
    __doc__ = GenModDoc.construct_module_doc('OutputName')

    _input_ports = [('run_name', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'[""]', 'optional':True}),
                                 ('subfolder_name', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'[""]', 'optional':True}),
                                 ('delete_previous', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["False"]', 'optional':True}), ]


    _output_ports = [('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)')]

    def compute(self):
        port_map = {'run_name': ('runname', None, True),
                    'subfolder_name': ('subfolder_name', None, True),
                    'delete_previous': ('delete_previous', None, True), }

        name_info = utils.map_ports(self, port_map)

        if not name_info.has_key('runname') and not name_info.has_key('subfolder_name'):
            raise ModuleError(self, "either 'run_name' or 'subfolder_name' must be supplied")

        if name_info['runname'] and not name_info['runname'].isalnum():
            raise ModuleError(self, "run_name cannot contain spaces or any characters other than letters and numbers")

        if name_info['delete_previous']:
            #  do our best to clear out any previous contents with this name
            if name_info['subfolder_name'] != "":
                subfolder = os.path.join(utils.getrootdir(), name_info['subfolder_name'])
                shutil.rmtree(subfolder, ignore_errors=True)

            if name_info['runname'] != "":
                for fname in os.listdir(utils.getrootdir()):
                    if "_" + name_info['runname'] + "_" in fname:
                        fname = os.path.join(utils.getrootdir(), name_info['runname'])
                        os.unlink(fname)
        subfolder = os.path.join(utils.getrootdir(), name_info['subfolder_name'])
        if name_info['subfolder_name'] != "" and not os.path.exists(subfolder):
            os.makedirs(subfolder)

        self.set_output('run_name_info', name_info)

class MDSBuilder(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('MDSBuilder')

    _input_ports = [('RastersWithPARCInfoCSV', '(gov.usgs.sahm:RastersWithPARCInfoCSV:Other)'),
                                 ('fieldData', '(gov.usgs.sahm:FieldData:DataInput)'),
                                 ('backgroundPointCount', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                                 # ('backgroundProbSurf', '(gov.usgs.sahm:BackgroundSurfaceGenerator:Tools)'),
                                 ('backgroundProbSurf', '(edu.utah.sci.vistrails.basic:File)'),
                                 ('Seed', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["{}"]'.format(utils.get_seed()), 'optional':True}),
                                 ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}), ]


    _output_ports = [('mdsFile', '(gov.usgs.sahm:MergedDataSet:Other)')]

    def compute(self):
        port_map = {'fieldData': ('fieldData', None, False),
                    'backgroundPointCount': ('pointCount', None, False),
                    'backgroundProbSurf': ('probSurfacefName', None, False),
                    'Seed': ('seed', utils.get_seed, True),
                    'run_name_info': ('run_name_info', None, False), }

        MDSParams = utils.map_ports(self, port_map)

        inputs_csvs = self.force_get_input_list('RastersWithPARCInfoCSV')
        if len(inputs_csvs) == 0:
            raise ModuleError(self, "Must supply at least one 'RastersWithPARCInfoCSV'/nThis is the output from the PARC module")
        if not type(inputs_csvs[0]) == str:
            inputs_csvs = [i.name for i in inputs_csvs]

        run_name_info = MDSParams.get('run_name_info')
        if run_name_info:
            if not type(run_name_info) == dict:
                run_name_info = run_name_info.contents
            subfolder = run_name_info.get('subfolder_name', "")
            runname = run_name_info.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(
                                                MDSParams.get('fieldData', ''))
            if subfolder == '' and runname == '':
                subfolder, runname = utils.get_previous_run_info(
                                        os.path.split(inputs_csvs[0])[0])

        key_inputs = []
        for input in ['fieldData']:
            if MDSParams.has_key(input):
                key_inputs.append(MDSParams[input])

        inputs_names = [utils.get_relative_path(f, self) for f in inputs_csvs]
        for fname in inputs_names:
            key_inputs.append(fname)

        if MDSParams.has_key('probSurfacefName'):
            key_inputs.append(MDSParams['probSurfacefName'])
        key_inputs.append(MDSParams['seed'])

        MDSParams['outputMDS'], signature, already_run = utils.make_next_file_complex(self,
                                        prefix='MergedDataset', suffix='.csv',
                                        key_inputs=key_inputs,
                                        subfolder=subfolder, runname=runname)

        if already_run:
            writetolog("No change in inputs or parameters using previous run of MDS Builder", True)
        else:
            inputs_csv = utils.mknextfile(prefix='CombinedPARCFiles', suffix='.csv', subfolder=subfolder, runname=runname)
            inputs_names = [utils.get_relative_path(f, self) for f in inputs_csvs]
            utils.merge_inputs_csvs(inputs_names, inputs_csv)
            MDSParams['inputsCSV'] = inputs_csv

            ourMDSBuilder = MDSB.MDSBuilder()
            utils.PySAHM_instance_params(ourMDSBuilder, MDSParams)

            writetolog("    inputsCSV=" + ourMDSBuilder.inputsCSV, False, False)
            writetolog("    fieldData=" + ourMDSBuilder.fieldData, False, False)
            writetolog("    outputMDS=" + ourMDSBuilder.outputMDS, False, False)

            try:
                ourMDSBuilder.run()
            except TrappedError as e:
                raise ModuleError(self, e.message)
            except:
                utils.informative_untrapped_error(self, "MDSBuilder")

        output_file = MDSParams['outputMDS']
        utils.write_hash_entry_pickle(signature, output_file)
        self.set_output('mdsFile', output_file)

class FieldDataQuery(SAHMDocumentedModule, Module):
    '''
    A wrapper to instantiate and run the FieldDataQuery module from PySAHM
    '''
    __doc__ = GenModDoc.construct_module_doc('FieldDataQuery')

    _input_ports = [('fieldData_file', '(gov.usgs.sahm:FieldData:DataInput)'),
                                 ('x_column', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["1"]', 'optional':True}),
                                 ('y_column', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["2"]', 'optional':True}),
                                 ('Response_column', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["3"]', 'optional':True}),
                                 ('Response_Presence_value', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["1"]', 'optional':True}),
                                 ('Response_Absence_value', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["0"]', 'optional':True}),
                                 ('ResponseType', '(edu.utah.sci.vistrails.basic:String)',
                                    {'entry_types': "['enum']",
                                        'values': "[['Presence(Absence)', 'Count']]", 'optional': True,
                                        'defaults':'["Presence(Absence)"]'}),
                                  ('Query_column', '(edu.utah.sci.vistrails.basic:String)', {'optional':True}),
                                  ('Query', '(edu.utah.sci.vistrails.basic:String)', {'optional':True}),
                                  ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}), ]
    _output_ports = [('fieldData', '(gov.usgs.sahm:FieldData:DataInput)'), ]

    def compute(self):
        writetolog("\nRunning FieldDataQuery", True)
        port_map = {'fieldData_file': ('fieldData', None, True),
            'x_column': ('x_col', None, True),
            'y_column': ('y_col', None, True),
            'Response_column': ('res_col', None, True),
            'Response_Presence_value': ('res_pres_val', None, True),
            'Response_Absence_value': ('res_abs_val', None, True),
            'ResponseType': ('response_type', None, True),
            'Query_column': ('query_col', None, False),
            'Query': ('query', None, False),
            'run_name_info': ('run_name_info', None, False), }

        FDQParams = utils.map_ports(self, port_map)
#          FDQOutput = utils.mknextfile(prefix='FDQ_', suffix='.csv')

        infile = open(FDQParams['fieldData'], "rb")
        csvReader = csv.DictReader(infile)

        run_name_info = FDQParams.get('run_name_info')
        if run_name_info:
            if not type(run_name_info) == dict:
                run_name_info = run_name_info.contents
            subfolder = run_name_info.get('subfolder_name', "")
            runname = run_name_info.get('runname', "")
        else:
            subfolder, runname = "", ""

        FDQOutput, signature, already_run = utils.make_next_file_complex(self,
                                        prefix='FDQ', suffix='.csv',
                                        key_inputs=[FDQParams['fieldData']],
                                        subfolder=subfolder, runname=runname)

        if already_run:
            writetolog("No change in inputs or parameters using previous run of FieldDataQuery", True)
        else:
            outfile = open(FDQOutput, "wb")
            csvwriter = csv.writer(outfile)
            if FDQParams["response_type"] == 'Count':
                responsetype = 'responseCount'
            else:
                responsetype = 'responseBinary'

            csvwriter.writerow(['X', 'Y', responsetype, "input=" + FDQParams['fieldData']])

            header = csvReader.fieldnames
            x_key = self.find_column(header, FDQParams['x_col'])
            y_key = self.find_column(header, FDQParams['y_col'])
            res_key = self.find_column(header, FDQParams['res_col'])

            use_query = False
            if self.has_input('Query'):
                use_query = True
                query = FDQParams['query']
                #  check if we're using a simple (equality) or complex (python syntax) query
                use_complex = any(s in query for s in ['[' + s + ']' for s in header])

            if self.has_input('Query_column'):
                query_col_key = self.find_column(header, FDQParams['query_col'])
            else:
                query_col_key = None

            for row in csvReader:
                if not use_query:
                    include_row = True
                elif use_complex:
                    include_row = self.complex_query(row, query)
                else:
                    include_row = self.simple_query(row, query, query_col_key)

                if include_row:
                    response = row[res_key]
                    if response.lower() in ["1", "true", "t", "present", "presence", str(FDQParams['res_pres_val']).lower()]:
                        response = 1
                    elif response.lower() in ["0", "false", "f", "absent", "absense", str(FDQParams['res_abs_val']).lower()]:
                        response = 0
                    elif responsetype == 'responseBinary':
                        try:
                            response = int(response)
                            if response > 0:
                                response = 1
                        except ValueError:
                            response = row[res_key]
                    else:
                        response = row[res_key]

                    csvwriter.writerow([row[x_key],
                                        row[y_key],
                                        response])

            del infile
            del outfile

        utils.write_hash_entry_pickle(signature, FDQOutput)
        self.set_output('fieldData', FDQOutput)


    def find_column(self, header, column):
        try:
            index = int(column) - 1
            if index > len(header) - 1:
                msg = "Field data input contains fewer columns than the number specified\n"
                msg += str(index + 1) + " is greater than " + str(len(header))
                writetolog(msg, True, True)
                raise ModuleError(self, msg)
            return header[index]
        except ValueError:
            if column in header:
                return column
            else:
                msg = "The specified column wasn't in the input file\n"
                msg += column + " not in " + str(header)
                writetolog(msg, True, True)
                raise ModuleError(self, msg)

    def simple_query(self, row, query, query_col):
        return row[query_col] == query

    def complex_query(self, row, query):

        for key in row.keys():
            query = query.replace('[' + key + ']', row[key])
        try:
            return eval(query)
        except NameError:
            msg = "There was a 'NameError' in the complex query you entered.\n"
            msg += "This is often an indication that strings are not being properly quoted in the python syntax.\n"
            msg += "Try enclosing the [fieldName] item in quotes.\n\n"
            msg += 'For example:  "[SourceType]" == "Expert"  instead of  [SourceType] == "Expert"'
            writetolog(msg, True, True)
            raise ModuleError(self, msg)

class FieldDataAggregateAndWeight(SAHMDocumentedModule, Module):
    '''
    Sanity!
    '''
    _input_ports = [('templateLayer', '(gov.usgs.sahm:TemplateLayer:DataInput)'),
                                 ('fieldData', '(gov.usgs.sahm:FieldData:DataInput)'),
                                 ('PointAggregationOrWeightMethod', '(edu.utah.sci.vistrails.basic:String)',
                                    {'entry_types': "['enum']",
                                     'values': "[['Collapse In Pixel', 'Weight Per Pixel']]", 'optional': True,
                                     'defaults':'["Collapse In Pixel"]'}),
                                 ('FD_EPSG_projection', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                                  ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}),
                                  ('drop_nodata_points', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}), ]
    _output_ports = [('fieldData', '(gov.usgs.sahm:FieldData:DataInput)')]

    __doc__ = GenModDoc.construct_module_doc('FieldDataAggregateAndWeight')

    def compute(self):
        writetolog("\nFieldDataAggregateAndWeight", True)
        port_map = {'templateLayer': ('templatefName', None, True),
            'fieldData': ('csv', None, True),
            'PointAggregationOrWeightMethod': ('aggMethod', None, True),
            'SDofGaussianKernel': ('sd', None, False),
            'FD_EPSG_projection': ('epsg', None, False),
            'run_name_info': ('run_name_info', None, False),
            'drop_nodata_points':('drop_nodata_points', None, True), }

        FDAWParams = utils.map_ports(self, port_map)
#          output_fname = utils.mknextfile(prefix='FDAW_', suffix='.csv')

        run_name_info = FDAWParams.get('run_name_info')
        if run_name_info:
            if not type(run_name_info) == dict:
                run_name_info = run_name_info.contents
            subfolder = run_name_info.get('subfolder_name', "")
            runname = run_name_info.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(FDAWParams['csv'])

        template_fname = FDAWParams['templatefName']
        if os.path.isdir(template_fname):
            template_fname = os.path.join(template_fname, "hdr.adf")
        output_fname, signature, already_run = utils.make_next_file_complex(self,
                                        prefix='FDAW', suffix='.csv',
                                        key_inputs=[FDAWParams['csv'], utils.get_raster_files(template_fname)],
                                        subfolder=subfolder, runname=runname)

        if already_run:
            writetolog("No change in inputs or parameters using previous run of FieldDataAggregateAndWeight", True)
        else:
            writetolog("    output_fname=" + output_fname, True, False)
            FDAWParams['output'] = output_fname

            ourFDAW = FDAW.FieldDataQuery()
            utils.PySAHM_instance_params(ourFDAW, FDAWParams)
            ourFDAW.processCSV()

        writetolog("Finished running FieldDataAggregateAndWeight", True)
        utils.write_hash_entry_pickle(signature, output_fname)
        self.set_output('fieldData', output_fname)

class PARC(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('PARC')

    _input_ports = [('predictor', "(gov.usgs.sahm:Predictor:DataInput)"),
                                ('PredictorList', '(gov.usgs.sahm:PredictorList:Other)'),
                                ('RastersWithPARCInfoCSV', '(gov.usgs.sahm:RastersWithPARCInfoCSV:Other)'),
                                ('templateLayer', '(gov.usgs.sahm:TemplateLayer:DataInput)'),
                                ('ignoreNonOverlap', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["False"]', 'optional':True}),
                                ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}), ]

    _output_ports = [('RastersWithPARCInfoCSV', '(gov.usgs.sahm:RastersWithPARCInfoCSV:Other)')]

    def compute(self):
        #  writetolog("\nRunning PARC", True)

        ourPARC = parc.PARC()
        template = utils.get_relative_path(self.force_get_input('templateLayer'), self)
        template_path, template_fname = os.path.split(template)
        template_fname = SpatialUtilities.getRasterShortName(template)

        run_name_info = self.force_get_input('run_name_info', None)
        if run_name_info:
            if not type(run_name_info) == dict:
                run_name_info = run_name_info.contents
            subfolder = run_name_info.get('subfolder_name', "")
            runname = run_name_info.get('runname', "")
            if runname:
                output_dname = os.path.join(utils.getrootdir(), subfolder, 'PARC_' + runname + "_" + template_fname)
            else:
                output_dname = os.path.join(utils.getrootdir(), subfolder, 'PARC_' + template_fname)
        else:
            subfolder, runname = "", ""
            output_dname = os.path.join(utils.getrootdir(), 'PARC_' + template_fname)

        if not os.path.exists(output_dname):
            os.mkdir(output_dname)

        if configuration.verbose:
            ourPARC.verbose = True
        ourPARC.logger = utils.getLogger()

        ourPARC.out_dir = output_dname

        ourPARC.processingMode = configuration.cur_processing_mode

        if self.has_input("ignoreNonOverlap"):
            ourPARC.ignoreNonOverlap = self.get_input("ignoreNonOverlap")

        key_inputs = utils.get_raster_files(template)
        for rasters_csv in self.force_get_input_list("RastersWithPARCInfoCSV"):
            key_inputs.append(utils.get_relative_path(rasters_csv))
        for predictor_list in self.force_get_input_list("PredictorList"):
            key_inputs.append(str(predictor_list))
        for predictor in self.force_get_input_list("predictor"):
            key_inputs.append(str(predictor))

        workingCSV, signature, already_run = utils.make_next_file_complex(self,
                                        prefix='PARCFiles', suffix='.csv',
                                        key_inputs=key_inputs,
                                        subfolder=os.path.join(subfolder, output_dname), runname=runname)


#          workingCSV = os.path.join(output_dname, "tmpFilesToPARC.csv")
        if already_run:
            writetolog("No change in inputs or parameters using previous run of PARC", True)
        else:
            f = open(workingCSV, "wb")
            csvWriter = csv.writer(f)
            csvWriter.writerow(["FilePath", "Categorical", "Resampling", "Aggregation"])

            if self.has_input("RastersWithPARCInfoCSV"):
                for input_rasters_csv in self.force_get_input_list('RastersWithPARCInfoCSV'):
                    input_rasters_csv_fname = utils.get_relative_path(input_rasters_csv)
                    csvReader = csv.reader(open(input_rasters_csv_fname), delimiter=",")
                    header = csvReader.next()
                    for row in csvReader:
                        csvWriter.writerow([utils.get_relative_path(row[0]), row[1], row[2], row[3]])


            if self.has_input("PredictorList"):
                predictor_lists = self.force_get_input_list('PredictorList')
                for predictor_list in predictor_lists:
                    for predictor in predictor_list:
                        csvWriter.writerow([utils.get_relative_path(predictor[0], self), predictor[1], predictor[2], predictor[3]])

            if self.has_input("predictor"):
                predictor_list = self.force_get_input_list('predictor')
                for predictor in predictor_list:
                    csvWriter.writerow([utils.get_relative_path(predictor[0], self), predictor[1], predictor[2], predictor[3]])
            f.close()
            del csvWriter
            ourPARC.inputs_CSV = workingCSV
            ourPARC.template = template
            writetolog('    template layer = ' + template)
            writetolog("    output_dname=" + output_dname, False, False)
            writetolog("    workingCSV=" + workingCSV, False, False)
            try:
                ourPARC.parcFiles()
            except TrappedError as e:
                writetolog(e.message)
                raise ModuleError(self, e.message)
            except:
                utils.informative_untrapped_error(self, "PARC")

        utils.write_hash_entry_pickle(signature, workingCSV)

        self.set_output('RastersWithPARCInfoCSV', workingCSV)

class Reclassifier(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('Reclassifier')

    _input_ports = [("inputRaster", "(edu.utah.sci.vistrails.basic:Path)"),
                    ('reclassFile', '(edu.utah.sci.vistrails.basic:File)', {'optional':True}),
                    ('reclassFileContents', '(edu.utah.sci.vistrails.basic:String)', {'optional':True}),
                    ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)')]

    _output_ports = [('outputRaster', '(edu.utah.sci.vistrails.basic:File)')]

    def compute(self):
        writetolog("\nRunning Reclassifier", True)
        port_map = {'inputRaster':('inputRaster', utils.dir_path_value, False),
                    'reclassFile':('reclassFile', utils.dir_path_value, False),
                    'reclassFileContents':('reclassFileContents', None, False), }

        argsDict = utils.map_ports(self, port_map)

        run_name_info = self.force_get_input('run_name_info', None)
        if run_name_info:
            if not type(run_name_info) == dict:
                run_name_info = run_name_info.contents
            subfolder = run_name_info.get('subfolder_name', "")
            runname = run_name_info.get('runname', "")
            if runname:
                output_dname = os.path.join(utils.getrootdir(), subfolder)
            else:
                output_dname = os.path.join(utils.getrootdir(), subfolder)
        else:
            subfolder, runname = "", ""
            output_dname = os.path.join(utils.getrootdir())

        from pySAHM.TiffProcessor import rasterReclassifier
        ourReclassifier = rasterReclassifier()
        ourReclassifier.inputFname = argsDict['inputRaster']

        if argsDict.has_key('reclassFileContents'):
            reclassFileName = utils.mknextfile("reclass", ".txt")
            reclassFile = open(reclassFileName, "w")
            reclassFile.write(self.force_get_input('reclassFileContents'))
            reclassFile.close()
            ourReclassifier.reclassFName = reclassFileName
        elif argsDict.has_key('reclassFile'):
            ourReclassifier.reclassFName = argsDict['reclassFile']
        else:
            msg = "Neither a reclass File or reclassFileContents have been specified\n"
            msg += "One or the other must be provided."
            raise ModuleError(self, msg)

        ourReclassifier.outDir = utils.getrootdir()

        in_shortname = utils.getShortName(ourReclassifier.inputFname)

        out_fname, signature, already_run = utils.make_next_file_complex(self,
                                        prefix=in_shortname, suffix='.tif',
                                        key_inputs=[argsDict['inputRaster']],
                                        subfolder=subfolder, runname=runname)

        if not already_run:
            ourReclassifier.outName = out_fname
            ourReclassifier.run()
            utils.write_hash_entry_pickle(signature, out_fname)

        self.set_output('outputRaster', out_fname)

class ReclassifierConfiguration(StandardModuleConfigurationWidget):
    #  FIXME add available_dict as parameter to allow config
    def __init__(self, module, controller, parent=None):

        StandardModuleConfigurationWidget.__init__(self, module, controller,
                                                   parent)
        self.setWindowTitle("Reclassification")
        self.build_gui()

        self.loadText()

    def build_gui(self):
        QtGui.QWidget.__init__(self)

        self.buttonSave = QtGui.QPushButton('Save', self)
        self.buttonReset = QtGui.QPushButton('Cancel', self)

        self.buttonSave.clicked.connect(self.handleSave)
        self.buttonReset.clicked.connect(self.handleReset)

        layout = QtGui.QVBoxLayout()
        self.textBox = QtGui.QTextEdit(self)

        layout.addWidget(self.textBox)

        buttonLayout = QtGui.QHBoxLayout()
        buttonLayout.addWidget(self.buttonSave)
        buttonLayout.addWidget(self.buttonReset)
        layout.addLayout(buttonLayout)
        self.setLayout(layout)

        self.path = None

    def getPortValue(self, portName):
        for i in xrange(self.module.getNumFunctions()):
            if self.module.functions[i].name == portName:
                return self.module.functions[i].params[0].strValue
        return None

    def handleSave(self):
        #  call this to save any current changes
        curStringValue = str(self.textBox.toPlainText())

        self.updateVisTrail(curStringValue)

    def handleReset(self):
        self.close()

#    def save(self):
#        with open(unicode(self.path), 'wb') as stream:
#            writer = csv.writer(stream)
#            #surely there is some cleaner way to get the header list!
#            header = [str(self.contents.horizontalHeaderItem(i).text())
#                    for i in range(self.contents.horizontalHeader().count())]
#            writer.writerow(header)
#            for row in range(self.contents.rowCount()):
#                rowdata = []
#                for column in range(self.contents.columnCount()):
#                    item = self.contents.item(row, column)
#                    if item is not None:
#                        rowdata.append(
#                            unicode(item.text()).encode('utf8'))
#                    else:
#                        rowdata.append('')
#                writer.writerow(rowdata)

    def updateVisTrail(self, strCurContents):
        self.controller.update_ports_and_functions(self.module.id,
                                           [], [], [("reclassFileContents", [strCurContents])])
        self.state_changed = False
        self.emit(QtCore.SIGNAL("stateChanged"))
        self.emit(QtCore.SIGNAL('doneConfigure'), self.module.id)

    def loadText(self):

        if self.getPortValue('reclassFileContents'):
            self.textBox.setText(self.getPortValue('reclassFileContents'))
        elif self.getPortValue('reclassFile'):
            curContents = open(self.getPortValue('reclassFile'), 'r').readlines()

class CategoricalToContinuous(Module):
    '''
    '''
#    __doc__ = GenModDoc.construct_module_doc('RasterFormatConverter')

    _input_ports = [("inputRaster", "(edu.utah.sci.vistrails.basic:File)"),
                    ('templateFile', '(gov.usgs.sahm:TemplateLayer:DataInput)'),
                    ]

    _output_ports = [('outputsPredictorListFile', '(gov.usgs.sahm:RastersWithPARCInfoCSV:Other)')]

    def compute(self):
        writetolog("\nCategoricalToContinuous", True)
        port_map = {'inputRaster':('inputRaster', utils.dir_path_value, True),
                    'templateFile':('templateFile', utils.dir_path_value, True)}

        argsDict = utils.map_ports(self, port_map)

        from pySAHM.TiffProcessor import categoricalToContinuousRasters
        ourC2C = categoricalToContinuousRasters()
        ourC2C.inputFname = argsDict['inputRaster']
        ourC2C.templateFName = argsDict['templateFile']
        shortName = SpatialUtilities.getRasterShortName(argsDict['inputRaster'])

        ourC2C.outDir = os.path.join(utils.getrootdir(), shortName + "_c2c")
        ourC2C.run()

        self.set_output('outputsPredictorListFile', ourC2C.outputPredictorsList)

class RasterFormatConverter(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('RasterFormatConverter')

    #  configuration = []
    _input_ports = [("inputMDS", "(gov.usgs.sahm:MergedDataSet:Other)"),
                    ('inputDir', '(edu.utah.sci.vistrails.basic:Directory)'),
                    ('format', '(edu.utah.sci.vistrails.basic:String)', {'optional':True}),
                    ('multipleCores', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True})]

    _output_ports = [('outputDir', '(edu.utah.sci.vistrails.basic:Directory)')]

    def compute(self):
        writetolog("\nRunning TiffConverter", True)
        ourRFC = RFC.FormatConverter()
        if self.has_input('inputMDS'):
            ourRFC.MDSFile = utils.get_relative_path(self.force_get_input('inputMDS'), self)
        elif self.has_input('inputDir'):
            ourRFC.inputDir = utils.get_relative_path(self.force_get_input('inputDir'), self)

        if self.has_input('format'):
            f = self.force_get_input('format')
            if f == '':
                f = 'asc'
            ourRFC.format = f

        if self.has_input("multipleCores"):
            if self.get_input("multipleCores"):
                ourRFC.multicores = "True"

        ourRFC.outputDir = utils.mknextdir(prefix='ConvertedRasters')
        if configuration.verbose:
            ourRFC.verbose = True
        ourRFC.logger = utils.getLogger()
        writetolog("    output directory = " + ourRFC.outputDir, False, False)

        try:
            ourRFC.run()
        except TrappedError as e:
            raise ModuleError(self, e.message)
        except:
            utils.informative_untrapped_error(self, "RasterFormatConverter")

        self.set_output('outputDir', ourRFC.outputDir)
        writetolog("\nFinished running TiffConverter", True)

class ModelEvaluationSplit(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('ModelEvaluationSplit')

    _input_ports = [("inputMDS", "(gov.usgs.sahm:MergedDataSet:Other)"),
                    ('trainingProportion', '(edu.utah.sci.vistrails.basic:Float)',
                        {'defaults':'["0.7"]', 'optional':True}),
                    ('Seed', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                    ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)')]
    _output_ports = [("outputMDS", "(gov.usgs.sahm:MergedDataSet:Other)")]

    def compute(self):
        writetolog("\nGenerating Model Evaluation split ", True)

        port_map = {'inputMDS':('i', utils.dir_path_value, True),
                'trainingProportion':('p', None, False),
                'Seed':('seed', None, False)}

        args = utils.map_ports(self, port_map)

        if self.has_input('run_name_info'):
            runinfo = self.force_get_input('run_name_info')
            if not type(runinfo) == dict:
                runinfo = runinfo.contents
            subfolder = runinfo.get('subfolder_name', "")
            runname = runinfo.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(args['i'])

        global models_path

        args['rc'] = utils.MDSresponseCol(args['i'])

        if args.get('p', 0.5) <= 0 or args.get('p', 0.5) > 1:
            raise ModuleError(self, "Train Proportion (trainProp) must be a number between 0 and 1 excluding 0")

        args['es'] = "TRUE"
        args['seed'] = utils.get_seed(args.get('seed', None))
        writetolog("    seed used for Split = " + str(args['seed']))

        outputMDS, signature, already_run = utils.make_next_file_complex(self,
                                prefix='ModelEvaluationSplit', suffix='.csv',
                                key_inputs=[args['i']],
                                subfolder=subfolder, runname=runname)
        args['o'] = outputMDS

        if not already_run:
            utils.run_R_script("TestTrainSplit.r", args, self, new_r_path=configuration.r_path)
            utils.write_hash_entry_pickle(signature, outputMDS)

        writetolog("Finished Model Evaluation split ", True)
        self.set_output("outputMDS", outputMDS)

class ModelSelectionSplit(SAHMDocumentedModule, Module):
    '''
    ToDo: Marian to write
    '''
    __doc__ = GenModDoc.construct_module_doc('ModelSelectionSplit')

    _input_ports = [("inputMDS", "(gov.usgs.sahm:MergedDataSet:Other)"),
                    ('trainingProportion', '(edu.utah.sci.vistrails.basic:Float)',
                        {'defaults':'["0.7"]', 'optional':True}),
                    ('Seed', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                    ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}), ]

    _output_ports = [("outputMDS", "(gov.usgs.sahm:MergedDataSet:Other)")]

    def compute(self):
        writetolog("\nGenerating Model Selection split ", True)
        port_map = {'inputMDS':('i', utils.dir_path_value, True),
                'trainingProportion':('p', None, False),
                'Seed':('seed', None, False)}

        args = utils.map_ports(self, port_map)

        if self.has_input('run_name_info'):
            runinfo = self.force_get_input('run_name_info')
            if not type(runinfo) == dict:
                runinfo = runinfo.contents
            subfolder = runinfo.get('subfolder_name', "")
            runname = runinfo.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(args['i'])

        global models_path
        args['rc'] = utils.MDSresponseCol(args['i'])

        if args.get('p', 0.5) <= 0 or args.get('p', 0.5) > 1:
            raise ModuleError(self, "Train Proportion (trainProp) must be a number between 0 and 1 excluding 0")

        args['es'] = "FALSE"
        args['seed'] = utils.get_seed(args.get('seed', None))
        writetolog("    seed used for Split = " + str(args['seed']))

        outputMDS, signature, already_run = utils.make_next_file_complex(self,
                                prefix='modelSelectionSplit', suffix='.csv',
                                key_inputs=[args['i']],
                                subfolder=subfolder, runname=runname)
        args['o'] = outputMDS

        if not already_run:
            utils.run_R_script("TestTrainSplit.r", args, self, new_r_path=configuration.r_path)

        writetolog("Finished Model Selection split ", True)
        utils.write_hash_entry_pickle(signature, outputMDS)
        self.set_output("outputMDS", outputMDS)

class ModelSelectionCrossValidation(SAHMDocumentedModule, Module):
    '''
    ToDo: Marian to write
    '''
    __doc__ = GenModDoc.construct_module_doc('ModelSelectionCrossValidation')

    _input_ports = [("inputMDS", "(gov.usgs.sahm:MergedDataSet:Other)"),
                    ('nFolds', '(edu.utah.sci.vistrails.basic:Integer)',
                        {'defaults':'["10"]', 'optional':True}),
                    ('SpatialSplit', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["False"]', 'optional':True}),
                    ('Stratify', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                    ('Seed', '(edu.utah.sci.vistrails.basic:Integer)', {'optional':True}),
                    ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}) ]
    _output_ports = [("outputMDS", "(gov.usgs.sahm:MergedDataSet:Other)")]

    def compute(self):
        writetolog("\nGenerating Cross Validation split ", True)
        port_map = {'inputMDS':('i', utils.dir_path_value, True),
                    'nFolds':('nf', None, True),
                    'SpatialSplit':('spt', utils.R_boolean, False),
                    'Stratify':('stra', utils.R_boolean, True)}

        argsDict = utils.map_ports(self, port_map)

        if self.has_input('run_name_info'):
            runinfo = self.force_get_input('run_name_info')
            if not type(runinfo) == dict:
                runinfo = runinfo.contents
            subfolder = runinfo.get('subfolder_name', "")
            runname = runinfo.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(argsDict['i'])

        argsDict["rc"] = utils.MDSresponseCol(argsDict["i"])

        if argsDict["nf"] <= 0:
            raise ModuleError(self, "Number of Folds must be greater than 0")

        seed = utils.get_seed(self.force_get_input("Seed", None))
        if not argsDict.has_key('spt'):
                argsDict['spt'] = 'FALSE'

        writetolog("    seed used for Split = " + str(seed))
        argsDict["seed"] = str(seed)
        outputMDS, signature, already_run = utils.make_next_file_complex(self,
                                prefix='modelSelectionCV', suffix='.csv',
                                key_inputs=[argsDict['i']],
                                subfolder=subfolder, runname=runname)
        argsDict["o"] = outputMDS

        if not already_run:
            utils.run_R_script("CrossValidationSplit.r", argsDict, self, new_r_path=configuration.r_path)

        writetolog("Finished Cross Validation split ", True)
        utils.write_hash_entry_pickle(signature, outputMDS)
        self.set_output("outputMDS", outputMDS)

class CovariateCorrelationAndSelection(SAHMDocumentedModule, Module):
    '''
    '''
    __doc__ = GenModDoc.construct_module_doc('CovariateCorrelationAndSelection')

    _input_ports = [("inputMDS", "(gov.usgs.sahm:MergedDataSet:Other)"),
                    ('selectionName', '(edu.utah.sci.vistrails.basic:String)', {'defaults':'["initial"]', 'optional':True}),
                    ('ShowGUI', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["True"]', 'optional':True}),
                    ('numPlots', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["8"]', 'optional':True}),
                    ('minCor', '(edu.utah.sci.vistrails.basic:Float)', {'defaults':'["0.7"]', 'optional':True}),
                    ('corsWithHighest', '(edu.utah.sci.vistrails.basic:Boolean)', {'defaults':'["False"]', 'optional':True}),
                    ('Seed', '(edu.utah.sci.vistrails.basic:Integer)', {'defaults':'["{}"]'.format(utils.get_seed()), 'optional':True}),
                    ('run_name_info', '(gov.usgs.sahm:OutputNameInfo:Other)', {'optional':False}), ]
    _output_ports = [("outputMDS", "(gov.usgs.sahm:MergedDataSet:Other)")]

    def compute(self):
        writetolog("\nOpening Select Predictors Layers widget", True)

        port_map = {'inputMDS': ('inputMDS', None, True),
                    'selectionName': ('selectionName', None, True),
                    'ShowGUI': ('ShowGUI', None, True),
                    'numPlots': ('numPlots', None, False),
                    'minCor': ('minCor', None, False),
                    'corsWithHighest': ('corsWithHighest', utils.R_boolean, False),
                    'Seed': ('seed', utils.get_seed, True)}

        params = utils.map_ports(self, port_map)

        if self.has_input('run_name_info'):
            runinfo = self.force_get_input('run_name_info')
            if not type(runinfo) == dict:
                runinfo = runinfo.contents
            subfolder = runinfo.get('subfolder_name', "")
            runname = runinfo.get('runname', "")
        else:
            subfolder, runname = utils.get_previous_run_info(params['inputMDS'])

        if runname:
            runname = runname + "_" + params['selectionName']
        else:
            runname = params['selectionName']

        writetolog("    seed used for subsampling = " + str(params['seed']))
        global session_dir

        outfname = os.path.join(session_dir, subfolder, "CovariateCorrelationOutputMDS_" + runname + ".csv")

        if outfname == params['inputMDS']:
            outfname = outfname[:-4] + "_2.csv"

        params['outputMDS'] = outfname
        params['displayJPEG'] = os.path.join(session_dir, subfolder, "CovariateCorrelationDisplay.png")
        params['r_path'] = configuration.r_path
        params['module'] = self
        writetolog("    inputMDS = " + params['inputMDS'], False, False)
        writetolog("    displayJPEG = " + params['displayJPEG'], False, False)
        writetolog("    outputMDS = " + params['outputMDS'], False, False)

        if os.path.exists(params['outputMDS']) and params['ShowGUI']:
            utils.applyMDS_selection(params['outputMDS'], params['inputMDS'])
            os.remove(params['outputMDS'])
            self.callDisplayMDS(params)
        elif os.path.exists(params['outputMDS']) and not params['ShowGUI']:
            utils.applyMDS_selection(params['outputMDS'], params['inputMDS'])
            os.remove(params['outputMDS'])
            shutil.copy2(params['inputMDS'], params['outputMDS'])
            writetolog("    Applying previous selection but not showing GUI", False, True)
        else:
            self.callDisplayMDS(params)


        output_file = params['outputMDS']
        writetolog("Finished Select Predictors Layers widget", True)
        self.set_output("outputMDS", output_file)

    def callDisplayMDS(self, kwargs):
        dialog = SelectListDialog(kwargs)
        #  dialog.setWindowFlags(QtCore.Qt.WindowMaximizeButtonHint)
        retVal = dialog.exec_()
        #  outputPredictorList = dialog.outputList
        if retVal == 1:
            raise ModuleError(self, "Cancel or Close selected (not OK) workflow halted.")

def load_max_ent_params():
    maxent_fname = os.path.join(os.path.dirname(__file__), 'maxent.csv')
    csv_reader = csv.reader(open(maxent_fname, 'rU'))
    #  pass on header
    csv_reader.next()
    input_ports = list(MAXENT._input_ports)

    docs = {}
    basic_pkg = 'edu.utah.sci.vistrails.basic'
    for row in csv_reader:
        [name, flag, p_type, default, doc, notes] = row
        name = name.strip()
        p_type = p_type.strip()
        kwargs = {}
        default = default.strip()
        if default:
            default = eval(default)
            kwargs['defaults'] = str([str(default)])
        kwargs['optional'] = True
        input_ports.append((name, '(' + basic_pkg + ':' + p_type + ')', kwargs))
        docs[name] = doc


    #  print 'MAXENT:', input_ports
    MAXENT._input_ports = input_ports
    MAXENT._port_docs = docs

    def provide_input_port_documentation(cls, port_name):
        return cls._port_docs[port_name]
    MAXENT.provide_input_port_documentation = \
        classmethod(provide_input_port_documentation)

def initialize():
    global maxent_path, java_path, color_breaks_csv
    global session_dir

    session_dir = configuration.cur_session_folder
    if not os.path.exists(session_dir):
        import tempfile
        orig_session_dir = session_dir
        session_dir = tempfile.mkdtemp(prefix="SAHM_session_dir_")
        utils.createLogger(session_dir, configuration.verbose)
        writetolog("!" * 79)
        writetolog("The previous session directory: " + orig_session_dir + " no longer exists on the file system!")
        writetolog("Defaulting to a random temporary location: " + session_dir)
        writetolog("!" * 79)

    utils.setrootdir(session_dir)
    utils.importOSGEO()
    utils.createLogger(session_dir, configuration.verbose)

    if system.systemType in ['Microsoft', 'Windows']:
        #  we're on a windows box, they probably installed VisTrails, R, and SAHM
        #  using the prebuilt msi,  let's make sure they're using the supplied R
        base_dir = os.path.split(os.path.split(sys.executable)[0])[0]
        default_r_dname = r"Central_R\R-3.5.0\bin"
        r_path = os.path.join(base_dir, default_r_dname)
        if os.path.exists(r_path):
            configuration.r_path = r_path
            configuration.set_deep_value('r_path', configuration.r_path)

            package_manager = get_package_manager()
            package = package_manager.get_package(identifier)
            package.persist_configuration()

    utils.set_r_path(configuration.r_path)

    try:
        testfname = os.path.join(utils.get_r_path(), "CanSAHMWriteToR.txt")
        open(testfname, "wb")
        os.remove(testfname)
    except:
        msg = ("!"*79 + "\n") * 3
        msg += "The current directory that R  is installed in:\n\t"
        msg += utils.get_r_path()
        msg += "\nIs not writeable!  This will cause errors in the\n"
        msg += "R modules unless all required packages are already installed!!!\n"
        msg += "Either point to an installation of R that is writeable or \n"
        msg += "Run VisTrails as administrator until all R packages have been downloaded.\n"
        msg += "\n  See page 3 of the user manual for more information!\n"
        msg += ("!"*79 + "\n") * 3
        writetolog(msg, True, True)


    maxent_path = os.path.abspath(configuration.maxent_path)
    if not os.path.exists(maxent_path) and os.path.exists(r"C:\Maxent\maxent.jar"):
        maxent_path = r"C:\Maxent"
        configuration.maxent_path = maxent_path
    if not os.path.exists(maxent_path) and maxent_path == r"..\\..\\Central_Maxent":
        maxent_path = r"C:\Maxent\maxent.jar"

    if not os.path.exists(maxent_path):
        msg = ("!"*79 + "\n") * 3
        msg += "The current installation of Maxent could not be found:\n\t"
        msg += maxent_path
        msg += "\nThe Maxent model will not work until this has been set correctly!\n"
        msg += "\n  See page 5 of the user manual for more information!\n"
        msg += ("!"*79 + "\n") * 3
        writetolog(msg, True, True)

    java_path = utils.find_java_exe(configuration.java_path)

    utils.set_seed(configuration.default_seed)

    gdal_data = os.path.join(os.path.dirname(__file__), "GDAL_Resources", "gdal-data")
    os.environ['GDAL_DATA'] = gdal_data
    projlib = os.path.join(os.path.dirname(__file__), "GDAL_Resources", "projlib")
    os.environ['PROJ_LIB'] = projlib

    color_breaks_csv = os.path.abspath(os.path.join(os.path.dirname(__file__), "ColorBreaks.csv"))

    load_max_ent_params()

    utilities.storeUNCDrives()
    utilities.start_new_pool(utilities.get_process_count(configuration.cur_processing_mode))

    global layers_csv_fname

    writetolog("*" * 79)
    writetolog("Initializing:", True, True)
    writetolog("  Locations of dependencies")
    writetolog("   Layers CSV = " + layers_csv_fname)
    writetolog("   R path = " + utils.get_r_path())
    writetolog("   Maxent folder = " + maxent_path)
    writetolog("    ")
    writetolog("*" * 79)
    writetolog("*" * 79)
    writetolog("SAHM output directory:   " + session_dir)
    writetolog("*" * 79)
    writetolog("*" * 79)

def finalize():
    pass

def generate_namespaces(modules):
    module_list = []
    for namespace, m_list in modules.iteritems():
        for module in m_list:
            m_dict = {'namespace': namespace}
            if type(module) == tuple:
                m_dict.update(module[1])
                module_list.append((module[0], m_dict))
                #  print 'm_dict:', m_dict
            else:
                module_list.append((module, m_dict))
    return module_list

def build_available_trees():
    trees = {}
    global layers_csv_fname
    layers_csv_fname = os.path.join(os.path.dirname(__file__), 'layers.csv')
    csv_reader = csv.reader(open(layers_csv_fname, 'rU'))
    csv_reader.next()
    first_file = csv_reader.next()[0]

    #  if the first file in the layers file does not exist assume that none
    #  of them do and use the exampledata version
    global atFORT

    atFORT = os.path.exists(first_file)

    if atFORT:
        csv_reader = csv.reader(open(layers_csv_fname, 'rU'))
        #  pass on header
        csv_reader.next()
        for row in csv_reader:
            if row[2] not in trees:
                trees[row[2]] = {}
            available_dict = trees[row[2]]
            if row[3] not in available_dict:
                available_dict[row[3]] = []
            available_dict[row[3]].append((row[0], row[1], row[4]))

    return trees

def build_predictor_modules():
    available_trees = build_available_trees()
    modules = []
    for name, tree in available_trees.iteritems():
        name_arr = name.strip().split()
        class_base = ''.join(n.capitalize() for n in name_arr)
        widget_class = get_predictor_widget(class_base, tree)
        config_class = get_predictor_config(class_base, tree)
        class_name = class_base + "Predictors"
        def get_widget_method(w_class):
            @staticmethod
            def get_widget_class():
                return w_class
            return get_widget_class
        module = type(class_name, (PredictorList,),
                      {'get_widget_class': get_widget_method(widget_class),
                       '_input_ports': \
                           [('value',
                             '(gov.usgs.sahm:%s:DataInput)' % class_name, True)]})
        modules.append((module, {'configureWidgetType': config_class,
                                 'moduleColor':INPUT_COLOR,
                                 'moduleFringe':INPUT_FRINGE}))
    for module in modules:
        module[0]._output_ports.append(('value_as_string', '(edu.utah.sci.vistrails.basic:String)', True))

    return modules

###################################

class TextFile(File):
    pass

class TextFileConfiguration(StandardModuleConfigurationWidget):
    #  FIXME add available_dict as parameter to allow config
    def __init__(self, module, controller, contents=None,
                 filter='', parent=None):
        StandardModuleConfigurationWidget.__init__(self, module, controller,
                                                   parent)
        self.fileFilter = filter

        if contents:
            self.contents = contents
        else:
            self.contents = QtGui.QTextEdit(self)

        self.setWindowTitle("Text File")
        self.build_gui()

        fid = self.findSourceFunction()
        if fid != -1:
            f = self.module.functions[fid]
            self.path = f.params[0].strValue
            self.loadText()

    def findSourceFunction(self):
        fid = -1
        for i in xrange(self.module.getNumFunctions()):
            if self.module.functions[i].name == "value":
                fid = i
                break
        return fid

    def build_gui(self):
        QtGui.QWidget.__init__(self)

        self.buttonOpen = QtGui.QPushButton('Open', self)
        self.buttonSave = QtGui.QPushButton('Save', self)
        self.buttonSaveAs = QtGui.QPushButton('Save As...', self)
        self.buttonReset = QtGui.QPushButton('Cancel', self)

        self.buttonOpen.clicked.connect(self.handleOpen)
        self.buttonSave.clicked.connect(self.handleSave)
        self.buttonSaveAs.clicked.connect(self.handleSaveAs)
        self.buttonReset.clicked.connect(self.handleReset)

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.contents)
        buttonLayout = QtGui.QHBoxLayout()
        buttonLayout.addWidget(self.buttonOpen)
        buttonLayout.addWidget(self.buttonSave)
        buttonLayout.addWidget(self.buttonSaveAs)
        buttonLayout.addWidget(self.buttonReset)
        layout.addLayout(buttonLayout)
        self.setLayout(layout)

        self.path = ''

    def handleReset(self):
        self.loadText()

    def handleSave(self):
        if not os.path.exists(self.path):
            self.path = QtGui.QFileDialog.getSaveFileName(
                self, 'Save File', os.path.split(self.path)[0], self.fileFilter)
        if not self.path.isEmpty():
            self.save()

    def save(self):
        f = open(self.path, "w")
        f.write(self.contents.toPlainText())

    def handleSaveAs(self):
        self.path = QtGui.QFileDialog.getSaveFileName(
                self, 'Save File As', os.path.split(self.path)[0], self.fileFilter)

        if not self.path.isEmpty():
            tmp = open(self.path, "w")
            del tmp
            self.handleSave()
            self.updateVisTrail()

    def updateVisTrail(self):
        self.controller.update_ports_and_functions(self.module.id,
                                           [], [], [("value", [str(self.path)])])
        self.state_changed = False
        self.emit(QtCore.SIGNAL("stateChanged"))
        self.emit(QtCore.SIGNAL('doneConfigure'), self.module.id)

    def handleOpen(self):
        self.path = QtGui.QFileDialog.getOpenFileName(
                self, 'Open File', os.path.split(self.path)[0], '')
        if not self.path.isEmpty():
            self.loadText()
            self.updateVisTrail()

    def loadText(self):
        f = open(self.path, 'r')
        data = f.read()
        self.contents.setText(data)

class CSVTextFile(TextFile):
    pass

class CSVTextFileConfiguration(TextFileConfiguration):
    #  FIXME add available_dict as parameter to allow config
    def __init__(self, module, controller, parent=None):

        fileFilter = 'CSV(*.csv)'
        contents = QtGui.QTableWidget(0, 0)
        TextFileConfiguration.__init__(self, module, controller, contents,
                                            fileFilter, parent)

        self.setWindowTitle("CSV Text File")


    def save(self):
        with open(unicode(self.path), 'wb') as stream:
            writer = csv.writer(stream)
            #  surely there is some cleaner way to get the header list!
            header = [str(self.contents.horizontalHeaderItem(i).text()) for i in
                      range(self.contents.horizontalHeader().count())]
            writer.writerow(header)
            for row in range(self.contents.rowCount()):
                rowdata = []
                for column in range(self.contents.columnCount()):
                    item = self.contents.item(row, column)
                    if item is not None:
                        rowdata.append(
                            unicode(item.text()).encode('utf8'))
                    else:
                        rowdata.append('')
                writer.writerow(rowdata)

    def loadText(self):
        with open(unicode(self.path), 'rb') as stream:
            csvReader = csv.reader(stream)
            header = csvReader.next()
            self.contents.setRowCount(0)
            self.contents.setColumnCount(len(header))
            self.contents.setHorizontalHeaderLabels(header)

            for rowdata in csvReader:
                row = self.contents.rowCount()
                self.contents.insertRow(row)
                self.contents.setColumnCount(len(rowdata))
                for column, data in enumerate(rowdata):
                    item = QtGui.QTableWidgetItem(data.decode('utf8'))
                    self.contents.setItem(row, column, item)

###################################

INPUT_COLOR = (0.76, 0.76, 0.8)
INPUT_FRINGE = [(0.0, 0.0),
                    (0.25, 0.0),
                    (0.0, 1.0)]

model_color = (0.76, 0.8, 0.76)
model_fringe = [(0.0, 0.0),
                    (0.25, 0.5),
                    (0.0, 1.0)]

output_color = (0.8, 0.8, 0.76)
output_fringe = [(0.0, 0.0),
                    (0.25, 0.0),
                    (0.0, 1.0)]

_modules = generate_namespaces({'DataInput': [
                                              (Predictor, {'moduleColor':INPUT_COLOR,
                                                           'moduleFringe':INPUT_FRINGE}),
                                              (PredictorListFile, {'moduleColor':INPUT_COLOR,
                                                           'moduleFringe':INPUT_FRINGE}),
                                              (FieldData, {'moduleColor':INPUT_COLOR,
                                                           'moduleFringe':INPUT_FRINGE}),
                                               (TemplateLayer, {'moduleColor':INPUT_COLOR,
                                                           'moduleFringe':INPUT_FRINGE}), ] + \
                                              build_predictor_modules(),
                                'Tools': [FieldDataQuery,
                                          FieldDataAggregateAndWeight,
                                          MDSBuilder,
                                          PARC,
                                          ModelEvaluationSplit,
                                          ModelSelectionSplit,
                                          ModelSelectionCrossValidation,
                                          CovariateCorrelationAndSelection,
                                          ApplyModel,
                                          BackgroundSurfaceGenerator,
                                          OutputName,
                                          EnsembleBuilder],
                                'GeospatialTools': [(Reclassifier, {'configureWidgetType': ReclassifierConfiguration}),
                                                    CategoricalToContinuous,
                                                    (GeoSpatialViewerCell, {'moduleColor':output_color,
                                                           'moduleFringe':output_fringe}),
                                                    (RasterLayer, {'moduleColor':INPUT_COLOR,
                                                           'moduleFringe':INPUT_FRINGE}),
                                                    (PolyLayer, {'moduleColor':INPUT_COLOR,
                                                           'moduleFringe':INPUT_FRINGE}),
                                                    (PointLayer, {'moduleColor':INPUT_COLOR,
                                                           'moduleFringe':INPUT_FRINGE}),
                                                    RasterFormatConverter,
#                                                      (LineLayer, {'moduleColor':INPUT_COLOR,
#                                                             'moduleFringe':INPUT_FRINGE}),
                                                    ],
                                'Models': [(GLM, {'moduleColor':model_color,
                                                           'moduleFringe':model_fringe}),
                                           (RandomForest, {'moduleColor':model_color,
                                                           'moduleFringe':model_fringe}),
                                           (MARS, {'moduleColor':model_color,
                                                           'moduleFringe':model_fringe}),
                                           (MAXENT, {'moduleColor':model_color,
                                                           'moduleFringe':model_fringe}),
                                           (BoostedRegressionTree,{'moduleColor':model_color,
                                                           'moduleFringe':model_fringe}),
                                           (UserDefinedCurve,{'moduleColor':model_color,
                                                           'moduleFringe':model_fringe}),
                                           (xgBoost,{'moduleColor':model_color,
                                                           'moduleFringe':model_fringe}),                                                           
                                           ],
                                'Other':  [(Model, {'abstract': True}),
                                           (VectorLayer, {'abstract': True}),
                                           (PredictorList, {'abstract': True}),
                                           (MergedDataSet, {'abstract': True}),
                                           (RastersWithPARCInfoCSV, {'abstract': True}),
                                           (BaseGeoViewerCell, {'abstract': True}),
                                           (OutputNameInfo, {'abstract': True}),
                                           (SAHMPathModule, {'abstract': True}),
                                           ],
                                'Output': [(ModelOutputViewer, {'moduleColor':output_color,
                                                           'moduleFringe':output_fringe}),
                                          (ModelMapViewer, {'moduleColor':output_color,
                                                           'moduleFringe':output_fringe}),
#                                            (ResponseCurveExplorer, {'moduleColor':output_color, # commented out until next release
#                                                             'moduleFringe':output_fringe}),
                                          ]
                                })

_upgrades = {}
_upgrades['DataInput|FieldData'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', None,
                                  dst_port_remap={'value': 'file'},
                                 src_port_remap={'value': 'value'})]
_upgrades['DataInput|TemplateLayer'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', None,
                                  dst_port_remap={'value': 'file'})]
_upgrades['DataInput|Predictor'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', None,
                                  dst_port_remap={'value': 'file'},
                                 function_remap={'AggregationMethod': utils.convert_old_enum,
                                              'ResampleMethod': utils.convert_old_enum, })]

_upgrades['Models|MAXENT'] = [UpgradeModuleRemap(None, '1.0.2', '1.0.2', None,
                                  dst_port_remap={'inputMDS': 'mdsFile'})]

for m in ['GLM', 'MARS', 'RandomForest', 'BoostedRegressionTree', 'MAXENT']:
        _upgrades['Models|' + m] = [UpgradeModuleRemap(None, '1.0.2', '1.0.2', 'Models|' + m,
                                        dst_port_remap={'modelWorkspace': utils.getParentDir}),
                                    UpgradeModuleRemap(None, '1.2.0', '1.2.0', 'Models|' + m,
                                        dst_port_remap={'UsePseudoAbs':None}),
                                    UpgradeModuleRemap(None, '2.0.0', '2.0.0', 'Models|' + m,
                                        dst_port_remap={'ThresholdOptimizationMethod': utils.convert_tom})]

_upgrades['Output|SAHMSpatialOutputViewerCell'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', 'Output|ModelMapViewer', {}), ]
_upgrades['Output|SAHMModelOutputViewerCell'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', 'Output|ModelOutputViewer', {})]

_upgrades['Tools|BackgroundSurfaceGenerator'] = [UpgradeModuleRemap(None, '1.0.2', '2.0.1', None,
                                  dst_port_remap={'bias': 'continuous'})]
_upgrades['Tools|MDSBuilder'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.1', None,
                                  dst_port_remap={'backgroundpointCount': 'backgroundPointCount',
                                                  'backgroundPointType':None})]
_upgrades['Tools|PARC'] = [UpgradeModuleRemap(None, '1.0.2', '1.0.2', None,
                                  dst_port_remap={'bias': '',
                                          'multipleCores': ''}),
                           UpgradeModuleRemap(None, '1.2.0', '2.0.0', None,
                                  dst_port_remap={'outputFolderName': None})]
_upgrades['Tools|RasterFormatConverter'] = [UpgradeModuleRemap(None, '1.0.2', '1.0.2', None,
                                  dst_port_remap={'multipleCores': ''})]
_upgrades['Tools|ApplyModel'] = [UpgradeModuleRemap(None, '1.0.1', '1.0.1', None,
                                 function_remap={'modelWorkspace': utils.getParentDir})]
_upgrades['Tools|FieldDataQuery'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', None,
                                  dst_port_remap={'ResponseType': 'ResponseType'},
                                 function_remap={'ResponseType': utils.convert_old_enum})]
_upgrades['Tools|FieldDataAggregateAndWeight'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', 'Tools|FieldDataAggregateAndWeight',
                      {'function_remap': {'ResponseType': utils.convert_old_enum} })]

_upgrades['GeospatialTools|RasterLayer'] = [UpgradeModuleRemap(None, '2.0.0', '2.0.0', None,
                                  dst_port_remap={'value': 'file'},
                                 function_remap={'cmap': utils.convert_old_enum})]
