import math
import os 
import numpy as np
import pandas as pd
import sqlite3

class SupportCalculations(object):
    
    conn = sqlite3.connect('tanks-4-09-data.db')
    
    def __init__(self, 
                 location, 
                 productMix, 
                 shellShade, 
                 shellCondition, 
                 roofShade, 
                 roofCondition, 
                 isShellInsulated, 
                 isRoofInsulated, 
                 shellHeight, 
                 tankDiameter):
        
        self.location = location
        self.productMix = productMix
        self.shellShade = shellShade
        self.shellCondition = shellCondition
        self.roofShade = roofShade
        self.roofCondition = roofCondition
        self.isShellInsulated = isShellInsulated
        self.isRoofInsulated = isRoofInsulated
        self.shellHeight = shellHeight
        self.tankDiameter = tankDiameter
        
    def resultsShellRoofData(self):
        
        queryStringShell1 = " FROM solar_abs_data WHERE Surface_Color == " + "'" + str(self.shellShade) + "'"
        queryStringShell = "SELECT " + str(self.shellCondition) + queryStringShell1
        sqlQueryqueryStringShell = pd.read_sql(queryStringShell, self.conn)
        alphaS_shell = sqlQueryqueryStringShell.iloc[0][self.shellCondition]

        queryStringRoof1 = " FROM solar_abs_data WHERE Surface_Color == " + "'" + str(self.roofShade) + "'"
        queryStringRoof = "SELECT " + str(self.roofCondition) + queryStringRoof1
        sqlQueryqueryStringRoof = pd.read_sql(queryStringRoof, self.conn)
        alphaR_roof = sqlQueryqueryStringRoof.iloc[0][self.roofCondition]
        
        return {'alphaS_shell': alphaS_shell, 'alphaR_roof': alphaR_roof}
        
    def SQLResultsLocation(self):
        
        locationString1 = self.location
        locationString2 = "SELECT * FROM revised_met_data WHERE location == "
        locationString3 = "'" + str(locationString1) + "'"
        queryStringLocation = locationString2 + locationString3
        sqlDataLocation = pd.read_sql(queryStringLocation, self.conn)
        
        return sqlDataLocation    
    
    def SQLResultsProductData(self):
                      
        _productMix = self.productMix
        
        if len(_productMix) == 1:
                
            if _productMix[0]['productClass'] in ['Crude Oils']:
                 #Does not accept mixture
                if _productMix[0]['productName'] in ['Crude oil (RVP 5)']:

                    _queryStringProduct1 = "SELECT * FROM chemical_data " 
                    _queryStringProduct2 = "WHERE CATEGORY == " + "'"  
                    _queryStringProduct3 = str(_productMix[0]['productClass']) + "'"
                    _queryStringProduct = _queryStringProduct1 + _queryStringProduct2 + _queryStringProduct3
                    sqlDataProduct = pd.read_sql(_queryStringProduct, self.conn)

                    return sqlDataProduct 

                else:
                    raise ValueError('Assumption Invalid')

            elif _productMix[0]['productClass'] in ['Petroleum Distillates']:
                 #Does not accept mixture
                if _productMix[0]['productName'] in ['Gasoline (RVP 10)', 'Gasoline (RVP 13)', 'Gasoline (RVP 7)',
                                              'Gasoline (RVP 6)', 'Gasoline (RVP 8)', 'Gasoline (RVP 9)',
                                              'Gasoline (RVP 11)', 'Gasoline (RVP 12)', 'Gasoline (RVP 7.8)',
                                              'Gasoline (RVP 8.3)', 'Gasoline (RVP 11.5)','Gasoline (RVP 13.5)']:

                    _queryStringProduct1 = "SELECT * FROM chemical_data " 
                    _queryStringProduct2 = "WHERE CATEGORY == " + "'" 
                    _queryStringProduct3 = str(_productMix[0]['productClass']) + "'"
                    _queryStringProduct4 = " AND NAME == " +  "'" + str(_productMix[0]['productName']) + "'"
                    _queryStringProduct = _queryStringProduct1 + _queryStringProduct2 + _queryStringProduct3 + _queryStringProduct4

                    sqlDataProduct = pd.read_sql(_queryStringProduct, self.conn)

                    return sqlDataProduct
                
                else:
                    raise ValueError('Assumption Invalid')
            
            else:
                raise ValueError('Assumption Invalid')
        
        elif _productMix[0]['productClass'] in ['Organic Liquids']:
            
            #Accepts Mixture
            
            nameList = [each['productName'] for each in _productMix]
            compList = [each['composition'] for each in _productMix]
            
            if sum(compList) >= 0.999999:

                strBuild = "("
                for idx, each in enumerate(nameList[:-1]):
                    strBuild += "'" + each + "'"  + ', '
                whereString = strBuild + "'" + nameList[-1]+ "'" + ")"

                queryStringProductOL1 = "SELECT * FROM chemical_data "
                queryStringProductOL2 = "WHERE name IN " + whereString
                queryStringProductOL3 = queryStringProductOL1 + queryStringProductOL2
                _sqlDataProduct = pd.read_sql(queryStringProductOL3, self.conn)

                sqlDataProduct = _sqlDataProduct.assign(composition = compList)

            return sqlDataProduct
        
        else:
            raise ValueError('Assumption Invalid Mixture')
            
    def locationData(self):
    
        locationData = self.SQLResultsLocation()
        
        alphaS_shell = self.resultsShellRoofData()['alphaS_shell']
        alphaR_roof = self.resultsShellRoofData()['alphaR_roof']
    
        tax_maxAmbientTemp_R = locationData[locationData['symbol'] == 'TAX'].values.tolist()[0][-1] + 460.67
        tan_minAmbientTemp_R = locationData[locationData['symbol'] == 'TAN'].values.tolist()[0][-1] + 460.67
        data_list_wind = locationData[locationData['symbol'] == 'V'].values.tolist()[0][-1]
        i_solarInsulation = locationData[locationData['symbol'] == 'I'].values.tolist()[0][-1]
        atmPressure = locationData[locationData['symbol'] == 'PA'].values.tolist()[0][-1]
        
        # Equation 1.11. Assumed to be the difference between the annual value.
        delta_ta_R = (tax_maxAmbientTemp_R - tan_minAmbientTemp_R)
        
        # Equation 1.30, value is in Rankine. Assumed to be the sum 
        #    of the annual values divided by 2
        taa_averageDailyTemp = (tan_minAmbientTemp_R + tax_maxAmbientTemp_R) * (1/2.)

        # Equation 1.31. Assumes that the tank is not insulated. 
        # Not sure about tb for fully insulated or partially insulated tanks.
        tb_liquidBulkTemp = taa_averageDailyTemp + (0.003*alphaS_shell*i_solarInsulation)
        
        return {'city': self.location,
                'shell_shade': self.shellShade,
                'shell_condition': self.shellCondition,
                'roof_shade': self.roofShade,
                'roof_condition': self.roofCondition,
                'isShellInsulated': self.isShellInsulated, 
                'isShellInsulated':self.isShellInsulated,
                'tax_maxAmbientTemp_R': tax_maxAmbientTemp_R, 
                'tan_minAmbientTemp_R': tan_minAmbientTemp_R, 
                'wind_speed': data_list_wind,
                'delta_ta_R':delta_ta_R, 
                'alphaS_shell': alphaS_shell, 
                'alphaR_roof': alphaR_roof, 
                'i_solarInsulation': i_solarInsulation,
                'taa_averageDailyTemp': taa_averageDailyTemp, 
                'atmPressure': atmPressure, 
                'tb_liquidBulkTemp': tb_liquidBulkTemp,                 
                'version': '06/2020, Table: 7.1-7'}
    
    def temperatureCalculations(self):
        
        locationTempData = self.locationData()
    
        #Fully Insulated
        if self.isShellInsulated == True and self.isRoofInsulated == True:
            
            # See Equation 1-5, note 1 for a fully insulated tank.
            # Assumes no cyclic heating of bulk liquid.
            deltatv_averageDailyVaporTempRange = 0
            
            # Equation 1-22, note 3, which references note 5.
            # TODO: Look into this issue for tb.
            # Issue here is that note 5 uses tb for a non-insulated tank, 
            # and we are working with a fully insulated tank.
            tla_averageDailyLiquidTemp = locationTempData['tb_liquidBulkTemp']
        
        #Partially Insulated
        elif self.isShellInsulated == True and self.isRoofInsulated == False:
            
            # Equation 1.8.
            # `delta_ta_R` is assumed to be in Rankine
            deltatv_averageDailyVaporTempRange = (0.6*locationTempData['delta_ta_R']) +\
                                                 (0.02*locationTempData['alphaR_roof'] *\
                                                  locationTempData['i_solarInsulation'])
            
            # Equation 1.29 for a partially insulated tank.
            # Again, the tb value is for a non-insulated tank.
            tla_averageDailyLiquidTemp = 0.3*locationTempData['taa_averageDailyTemp'] +\
                                            0.7*locationTempData['tb_liquidBulkTemp'] +\
                                            0.005*locationTempData['alphaR_roof'] *\
                                            locationTempData['i_solarInsulation']
        
        #Uninsulated
        elif self.isShellInsulated == False and self.isRoofInsulated == False:
            
            # Equation 1.6 for uninsulated tank.
            hs_d = self.shellHeight / self.tankDiameter
            
            del_tv_1 = (2.2 * hs_d) + 1.9
            del_tv_2 = 0.042*locationTempData['alphaR_roof']*locationTempData['i_solarInsulation']
            del_tv_3 = 0.026*hs_d*locationTempData['alphaS_shell']*locationTempData['i_solarInsulation']
            del_tv_4 = locationTempData['delta_ta_R']
            del_tv_5 = (1 - ((0.8)/(del_tv_1)))*del_tv_4
            del_tv_6 = (del_tv_2 + del_tv_3) / del_tv_1
            
            deltatv_averageDailyVaporTempRange = del_tv_5 + del_tv_6
            
            # Equation 1.27 for uninsulated tank.
            tla_1 = (4.4 * hs_d) + 3.8
            tla_2_1 = 0.021*locationTempData['alphaR_roof']*locationTempData['i_solarInsulation']
            tla_2_2 = 0.013*hs_d*locationTempData['alphaS_shell']*locationTempData['i_solarInsulation']
            
            pt_1 = (0.5-(0.8/tla_1))*locationTempData['taa_averageDailyTemp']
            pt_2 = (0.5+(0.8/tla_1))*locationTempData['tb_liquidBulkTemp']
            pt_3 = (tla_2_1 + tla_2_2) / tla_1
            
            tla_averageDailyLiquidTemp = pt_1 + pt_2 + pt_3
            
        else:
            raise ValueError('Assumption Invalid')
            
        # Figure 7.1-17, values are in Rankine
        tlx_averageDailyMaxLiqSurfaceTemp_R = tla_averageDailyLiquidTemp +\
                                            (0.25 * deltatv_averageDailyVaporTempRange)
        
        tln_averageDailyMinLiqSurfaceTemp_R = tla_averageDailyLiquidTemp -\
                                            (0.25 * deltatv_averageDailyVaporTempRange)
        
        return {'tlx_R': tlx_averageDailyMaxLiqSurfaceTemp_R,
                'tln_R': tln_averageDailyMinLiqSurfaceTemp_R, 
                'tla_R': tla_averageDailyLiquidTemp, 
                'delta_tv': deltatv_averageDailyVaporTempRange}
    
    def organicLiquidMixtureVaporMolecularWeight(self):
        
        _productMix = self.productMix
        
        if _productMix[0]['productClass'] in ['Organic Liquids']:
        
            productDataOL = self.SQLResultsProductData()
            tempData = self.temperatureCalculations()

            check_A = productDataOL['VP_COEF_A'].values.tolist()
            check_B = productDataOL['VP_COEF_B'].values.tolist()
            check_C = productDataOL['VP_COEF_C'].values.tolist()

            if 0 not in check_A:
                if 0 not in check_B:
                    if 0 not in check_C:

                        # tlx_averageDailyMaxLiqSurfaceTemp_R
                        productDataOL['log_pva'] = productDataOL['VP_COEF_A'] - ((productDataOL['VP_COEF_B']) /\
                                                                                       (((tempData['tla_R'] - 491.7)/1.8) + productDataOL['VP_COEF_C']))
                        
                        productDataOL['log_pvx'] = productDataOL['VP_COEF_A'] - ((productDataOL['VP_COEF_B']) /\
                                                                                       (((tempData['tlx_R'] - 491.7)/1.8) + productDataOL['VP_COEF_C']))
                        
                        productDataOL['log_pvn'] = productDataOL['VP_COEF_A'] - ((productDataOL['VP_COEF_B']) /\
                                                                                       (((tempData['tln_R'] - 491.7)/1.8) + productDataOL['VP_COEF_C']))
                        
                        productDataOL['pva'] = 10**(productDataOL['log_pva']) / 51.7149
                        productDataOL['pvx'] = 10**(productDataOL['log_pvx']) / 51.7149
                        productDataOL['pvn'] = 10**(productDataOL['log_pvn']) / 51.7149
                        
                        productDataOL['moles'] = 100*productDataOL['composition'] / productDataOL['MOLWT']
                        productDataOL['mol_frac'] =  productDataOL['moles'] / sum(productDataOL['moles'].values.tolist())
                        productDataOL['frac_pva'] = productDataOL['pva'] * productDataOL['composition']
                        productDataOL['frac_pvx'] = productDataOL['pvx'] * productDataOL['composition']
                        productDataOL['frac_pvn'] = productDataOL['pvn'] * productDataOL['composition']
                        
                        productDataOL['partial_pressure_pva'] = productDataOL['mol_frac'] * productDataOL['pva']
                        productDataOL['partial_pressure_pvx'] = productDataOL['mol_frac'] * productDataOL['pvx']
                        productDataOL['partial_pressure_pvn'] = productDataOL['mol_frac'] * productDataOL['pvn']
                        
                        total_pva = sum(productDataOL['partial_pressure_pva'].values.tolist())
                        total_pvx = sum(productDataOL['partial_pressure_pvx'].values.tolist())
                        total_pvn = sum(productDataOL['partial_pressure_pvn'].values.tolist())
                        
                        productDataOL['y_i_pva'] = productDataOL['partial_pressure_pva'] / total_pva
                        productDataOL['y_i_pvx'] = productDataOL['partial_pressure_pvx'] / total_pvx
                        productDataOL['y_i_pvn'] = productDataOL['partial_pressure_pvn'] / total_pvn
                        
                        productDataOL['mv_i_pva'] = productDataOL['MOLWT'] * productDataOL['y_i_pva'] 
                        productDataOL['mv_i_pvx'] = productDataOL['MOLWT'] * productDataOL['y_i_pvx'] 
                        productDataOL['mv_i_pvn'] = productDataOL['MOLWT'] * productDataOL['y_i_pvn'] 
        
                        vapor_mw = sum(productDataOL['mv_i_pva'].values.tolist())
            
                    else:
                        raise ValueError('Value C Not Present in Data')
                else:
                    raise ValueError('Value B Not Present in Data')
            else: 
                raise ValueError('Value A Not Present in Data')

        else:
            raise ValueError('Grant')
        
        return {'pva': total_pva, 
                'plx': total_pvx, 
                'pln': total_pvn, 
                'delta_pv': total_pvx - total_pvn,
                'vapor_mw': vapor_mw, 
                'rvp': None}
    
    def crudePetVaporMolecularWeight(self):
        
        _productMix = self.productMix
        
        if _productMix[0]['productClass'] in ['Crude Oils', 'Petroleum Distillates']:
            
            
            _temperatureCalculations = self.temperatureCalculations()
            _SQLResultsProductData = self.SQLResultsProductData()

            rvp = _SQLResultsProductData['REID']

            const_a = 12.82 - 0.9672*math.log(rvp)
            const_b = 7261 - 1216*math.log(rvp)

            constDict =  {'const_a': const_a, 'const_b': const_b}
            val_pva = constDict['const_a'] - ( (constDict['const_b']) / _temperatureCalculations['tla_R'] )
            val_pvx = constDict['const_a'] - ( (constDict['const_b']) / _temperatureCalculations['tlx_R'] )
            val_pvn = constDict['const_a'] - ( (constDict['const_b']) / _temperatureCalculations['tln_R'] )

            pva_trueVaporPressure = math.exp(val_pva)
            pvx_vapPressAveDailyMaxSurfaceTemp = math.exp(val_pvx)
            pvn_vapPressAveDailyMinSurfaceTemp = math.exp(val_pvn)
            vapor_mw = _SQLResultsProductData['VP_MOLWT'].values.tolist()[0]
            rvp = _SQLResultsProductData['REID'].values.tolist()[0]

            return {'pva': pva_trueVaporPressure, 
                    'plx': pvx_vapPressAveDailyMaxSurfaceTemp, 
                    'pln': pvn_vapPressAveDailyMinSurfaceTemp, 
                    'delta_pv': pvx_vapPressAveDailyMaxSurfaceTemp - pvn_vapPressAveDailyMinSurfaceTemp,
                    'vapor_mw': vapor_mw, 
                    'rvp': rvp}
        
        else:
            raise ValueError('Incorrect Assumption')
               
class FixedRoofTank(SupportCalculations):
    
    conn = sqlite3.connect('tanks-4-09-data.db')
    
    def __init__(self, 
                 location, 
                 productMix, 
                 shellShade, 
                 shellCondition, 
                 roofShade, 
                 roofCondition, 
                 isShellInsulated, 
                 isRoofInsulated, 
                 shellHeight, 
                 tankDiameter,
                 tankOrientation, 
                 isAverageLiquidHeightKnown, 
                 tankLength, 
                 roofType, 
                 tankDomeRoofRadius, 
                 roofSlope, 
                 breatherPressureSetting, 
                 breatherVaccumSetting):
        
        super().__init__(
                 location, 
                 productMix, 
                 shellShade, 
                 shellCondition, 
                 roofShade, 
                 roofCondition, 
                 isShellInsulated, 
                 isRoofInsulated, 
                 shellHeight, 
                 tankDiameter)
        
        self.tankOrientation = tankOrientation
        self.isAverageLiquidHeightKnown = isAverageLiquidHeightKnown
        self.tankLength = tankLength
        self.roofType = roofType 
        self.tankDomeRoofRadius = tankDomeRoofRadius
        self.roofSlope = roofSlope
        self.breatherPressureSetting = breatherPressureSetting
        self.breatherVaccumSetting = breatherVaccumSetting
        
    
    def stockVaporDensity(self):
        
        _locationData = self.locationData()
        _productMix = self.productMix
        hs_d = self.shellHeight / self.tankDiameter
        
        #Fully Insulated
        if self.isShellInsulated == True and self.isRoofInsulated == True:
            
            tv_averageVaporTemp = _locationData['tb_liquidBulkTemp']
            
        #Partially Insulated
        elif self.isShellInsulated == True and self.isRoofInsulated == False:
            
            tv_averageVaporTemp = 0.6 * _locationData['taa_averageDailyTemp'] + \
                                    0.4*_locationData['tb_liquidBulkTemp'] + \
                                    0.01 * _locationData['alphaR_roof'] * _loctionData['i_solarInsulation']
            
        #Partially Insulated
        elif self.isShellInsulated == False and self.isRoofInsulated == False:
            wv_pt_1 = (2.2 * hs_d + 1.1) * _locationData['taa_averageDailyTemp']
            wv_pt_2 = 0.8 * _locationData['tb_liquidBulkTemp']
            wv_pt_3 = 0.021 * _locationData['alphaR_roof'] * _locationData['i_solarInsulation']
            wv_pt_4 = 0.013 * hs_d * _locationData['alphaS_shell'] * _locationData['i_solarInsulation']
            wv_pt_5 = 2.2 * hs_d + 1.9
            
            tv_averageVaporTemp = (wv_pt_1 + wv_pt_2 + wv_pt_3 + wv_pt_4) / wv_pt_5
            
        else:
            raise ValueError('Error')
            
        if _productMix[0]['productClass'] in ['Crude Oils', 'Petroleum Distillates']:
            pva = self.crudePetVaporMolecularWeight()['pva']
            mw = self.crudePetVaporMolecularWeight()['vapor_mw']
            vaporDensity = (pva * mw) / (10.731 * tv_averageVaporTemp)
        
        elif _productMix[0]['productClass'] in ['Organic Liquids']:
            
            pva = self.organicLiquidMixtureVaporMolecularWeight()['pva']
            mw = self.organicLiquidMixtureVaporMolecularWeight()['vapor_mw']
            vaporDensity = (pva * mw) / (10.731 * tv_averageVaporTemp)
        
        else:
            raise ValueError('Grant')
            
        return vaporDensity
    
    def vaporSpaceOutage(self):
        '''
        Done - Check
        Equation 1-16. 
        Uses Notes 1, 2 and Equations 1-17 - 1-20.
        '''

        # Assumptions made within the function.
        tankShellRadius = self.tankDiameter / 2.
        
        if tankShellRadius >= 0.8*self.tankDiameter or\
            tankShellRadius <= 1.2*self.tankDiameter:
            pass
    
        # Assumed to be half of the shell height.
        if self.isAverageLiquidHeightKnown == True:
            aveLiquidLevel = self.averageLiquidHeight
            
        elif self.isAverageLiquidHeightKnown == False:
            aveLiquidLevel = self.shellHeight / 2.
            
        else:
            raise ValueError('Invalid Assumption')
            
        # hro_roof_outage
        if self.roofType == 'Dome':

            if self.isTankDomeRoofRadiusKnown == True:
                
                valueSqRoot = (self.tankDomeRoofRadius**2) -\
                              (tankShellRadius**2)
                
                if valueSqRoot >= 0:
                    #Equation 1.20.
                    roofHeight   =   self.tankDomeRoofRadius -\
                                        (valueSqRoot)**(0.5)
                
                else:
                    raise ValueError('Taking the negative of a square root')
                
                # Equation 1.19
                hro_roofOutage  =   roofHeight *\
                                    ( (1 / 2.) + ( (1 / 6.) *\
                                    (roofHeight / tankShellRadius)**2) )
            
            elif self.isTankDomeRoofRadiusKnown == False:
                
                #Why does this matter if the value is unknown?
                roofHeight = 0.268 * tankShellRadius 
                
                hro_roofOutage = 0.137 * tankShellRadius
            
            else:
                raise ValueError('Assumption Invalid')

        elif self.roofType == 'Cone':
            # Note 1 to determine hro_roof_outage for a cone roof.
            # Equation 1.17, and 1.18.
            hro_roofOutage = ( 1 / 3. ) * self.roofSlope * tankShellRadius
        
        else:
            raise ValueError('Incorrect Roof Type.')
        
        # hvo_vaporSpaceOutage
        if self.tankOrientation == 'Vertical':
            #Equation 1.16.
            hvo_vaporSpaceOutage =    self.shellHeight - \
                                        aveLiquidLevel + hro_roofOutage
            
            result = {'value': hvo_vaporSpaceOutage, 
                      'quantity': 'Vapor Space Outage (hvo)',
                      'equation': '1-16',
                      'version': '06/2020', 'elements': [{'hs': self.shellHeight, 
                                                          'hl': aveLiquidLevel, 
                                                          'hro': hro_roofOutage, 
                                                          'hr': self.roofSlope * tankShellRadius, 
                                                          'roof slope': self.roofSlope, 
                                                          'roof height': self.roofSlope * tankShellRadius,
                                                          'shell radius': tankShellRadius,
                                                          'shell height': self.shellHeight,
                                                          'liquid height': aveLiquidLevel,
                                                          'roof type': self.roofType,
                                                          'constant': None,
                                                          'd': self.tankDiameter,
                                                          'tank_orientation': self.tankOrientation,
                                                          'status': 'Done'}]} 
        elif self.tankOrientation == 'Horizontal':
            # See HVO, in Equation 1-16.
            hvo_vaporSpaceOutage = ( (math.pi) / 8. ) * self.tankDiameter
            
            result = {'value': hvo_vaporSpaceOutage, 
                      'quantity': 'Vapor Space Outage (hvo)',
                      'equation': '1-16',
                      'version': '06/2020', 
                      'elements': [{'hs': None, 
                                  'hl': None, 
                                  'hro': None, 
                                  'hr': None,
                                  'roof slope': None, 
                                  'roof height': None,
                                  'shell radius': None,
                                  'shell height': None,
                                  'liquid height': None,
                                  'roof type': None,
                                  'constant': 'pi/8',
                                  'd': self.tankDiameter,
                                  'tank_orientation': self.tankOrientation,
                                  'status': 'Done'}]} 
        
        else:
            raise ValueError('Incorrect Tank Orientation. \
                              Tank can either be Vertical \
                              or Horizontal.')
        
        return result
    
    def vaporSpaceVolume(self):
        '''
        Equation 1.3.
        '''
        if self.tankOrientation == 'Horizontal':

            effectiveDiameter = ( (self.tankLength * self.tankDiameter) /\
                                ( (1. / 4.) * (math.pi) ) ) ** (0.5)
            _tankDiameter = effectiveDiameter

        elif self.tankOrientation == 'Vertical':
            _tankDiameter = self.tankDiameter

        else:
            raise ValueError('Incorrect Tank Type')

        vv_vaporSpaceVolume =   (math.pi) * (1. / 4.) *\
                                (_tankDiameter ** 2) *\
                                self.vaporSpaceOutage()['value']

        return {'quantity': 'Tank Vapor Space Volume (vv)', 
                'equation': '1-3', 
                'value': vv_vaporSpaceVolume,
                'version': '06/2020',
                'notes': 'Uses effective diameter (eqn. 1-14) for horizontal tanks.',
                'elements': [{'constant': 'pi/4',
                              'tank orientation': self.tankOrientation,
                              'd': _tankDiameter, 
                              'hvo': self.vaporSpaceOutage()['value'], 
                              'status': 'Done'}]}
    
    def ventedVaporSpaceSatFactor(self):
        '''
        Partially Done - Outstanding item is to 
        compute pva (eqn 1-22) for mixtures.
        Equation 1-21.
        '''

        # Need to pass in a dict of all of the 
        # liquid components in the mixture, or a 
        # vapor pressure from a look up table.

        # TODO: Compute molecular weight of the vapor phase
        # TODO: Compuete total vapor pressure of the liquid.

        # Equation 1.22. See notes 1, and 2.
        # TODO: Build this part out for mixtures. What is here is for a single component.
        
        
        _productMix = self.productMix
        
        if _productMix[0]['productClass'] in ['Crude Oils', 'Petroleum Distillates']:
            vaporPresureArray = self.crudePetVaporMolecularWeight()
        
        elif _productMix[0]['productClass'] in ['Organic Liquids']:
            vaporPresureArray = self.organicLiquidMixtureVaporMolecularWeight()
        else:
            raise ValueError('Injcorrect Assumption')
        
        pva_vaporPressureAverageDaily = vaporPresureArray['pva']   

        # Equation 1.16
        hvo_vaporSpaceOutage = self.vaporSpaceOutage()['value']             

        ks_ventedVaporSpaceSatFactor = 1 / (1 + ( 0.053 *\
                                                 pva_vaporPressureAverageDaily *\
                                                 self.vaporSpaceOutage()['value'] ) )

        return {'quantity': 'Vented Vapor Saturation Factor (ks)', 
                'equation': '1-12', 
                'value': ks_ventedVaporSpaceSatFactor, 
                'version': '06/2020', 
                'notes': 'none', 
                'elements': [{'constant': 0.053, 
                              'pva': pva_vaporPressureAverageDaily, 
                              'hvo': hvo_vaporSpaceOutage, 
                              'status': 'Partially Done, need to include calculations for mixtures.'}]}
    
    def vaporSpaceExpansionFactor(self):
        '''
        Partially done -- Look into the value for tb 
        and the tank insulation, and the partial pressures. 
        Equation 1-5.
        
        Most of the values match the OK printout.
        
        Also, look into `delta_pv`, and need to fold in 
            logifc to pull in rvp based in user input.
    
        '''
        
        locationTempData = self.locationData()
        _productMix = self.productMix
        
        #Fully Insulated
        if self.isShellInsulated == True and self.isRoofInsulated == True:
            
            # See Equation 1-5, note 1 for a fully insulated tank.
            # Assumes no cyclic heating of bulk liquid.
            deltatv_averageDailyVaporTempRange = 0
            
            # Equation 1-22, note 3, which references note 5.
            # TODO: Look into this issue for tb.
            # Issue here is that note 5 uses tb for a non-insulated tank, 
            # and we are working with a fully insulated tank.
            tla_averageDailyLiquidTemp = locationTempData['tb_liquidBulkTemp']
        
        #Partially Insulated
        elif self.isShellInsulated == True and self.isRoofInsulated == False:
            
            # Equation 1.8.
            # `delta_ta_R` is assumed to be in Rankine
            deltatv_averageDailyVaporTempRange = (0.6*locationTempData['delta_ta_R']) +\
                                                 (0.02*locationTempData['alphaR_roof'] *\
                                                  locationTempData['i_solarInsulation'])
            
            # Equation 1.29 for a partially insulated tank.
            # Again, the tb value is for a non-insulated tank.
            tla_averageDailyLiquidTemp = 0.3*locationTempData['taa_averageDailyTemp'] +\
                                            0.7*locationTempData['tb_liquidBulkTemp'] +\
                                            0.005*locationTempData['alphaR_roof'] *\
                                            locationTempData['i_solarInsulation']
        
        #Uninsulated
        elif self.isShellInsulated == False and self.isRoofInsulated == False:
            
            # Equation 1.6 for uninsulated tank.
            hs_d = self.shellHeight / self.tankDiameter
            
            delta_tv_1 = (2.2 * hs_d) + 1.9
            delta_tv_2 = 0.042*locationTempData['alphaR_roof']*locationTempData['i_solarInsulation']
            delta_tv_3 = 0.026*hs_d*locationTempData['alphaS_shell']*locationTempData['i_solarInsulation']
            delta_tv_4 = locationTempData['delta_ta_R']
            
            delta_tv_5 = (1 - ((0.8)/(delta_tv_1)))* delta_tv_4
            delta_tv_6 = (delta_tv_2 + delta_tv_3) / delta_tv_1
            
            deltatv_averageDailyVaporTempRange = delta_tv_5 + delta_tv_6
            
            # Equation 1.27 for uninsulated tank.
            tla_1 = (4.4 * hs_d) + 3.8
            tla_2_1 = 0.021*locationTempData['alphaR_roof']*locationTempData['i_solarInsulation']
            tla_2_2 = 0.013*hs_d*locationTempData['alphaS_shell']*locationTempData['i_solarInsulation']
            
            pt_1 = (0.5-(0.8/tla_1))*locationTempData['taa_averageDailyTemp']
            pt_2 = (0.5+(0.8/tla_1))*locationTempData['tb_liquidBulkTemp']
            pt_3 = (tla_2_1 + tla_2_2) / tla_1
            
            tla_averageDailyLiquidTemp = pt_1 + pt_2 + pt_3
            
        else:
            raise ValueError('Assumption Invalid')
            
        # Figure 7.1-17, values are in Rankine
        tlx_averageDailyMaxLiqSurfaceTemp_R = tla_averageDailyLiquidTemp +\
                                            (0.25 * deltatv_averageDailyVaporTempRange)
        
        tln_averageDailyMinLiqSurfaceTemp_R = tla_averageDailyLiquidTemp -\
                                            (0.25 * deltatv_averageDailyVaporTempRange)
                
        # Only for Crude Oils, and Selected Petroleum Stocks
        if _productMix[0]['productClass'] in ['Crude Oils', 'Petroleum Distillates']:
            
            pva_trueVaporPressure = self.crudePetVaporMolecularWeight()['pva']
            pvx_vapPressAveDailyMaxSurfaceTemp = self.crudePetVaporMolecularWeight()['plx'] 
            pvn_vapPressAveDailyMinSurfaceTemp = self.crudePetVaporMolecularWeight()['pln'] 
        
        # TODO: Need to resolve this function for organic liquids.
        elif _productMix[0]['productClass'] in ['Organic Liquids']:
            
            pva_trueVaporPressure = self.organicLiquidMixtureVaporMolecularWeight()['pva']
            pvx_vapPressAveDailyMaxSurfaceTemp = self.organicLiquidMixtureVaporMolecularWeight()['plx'] 
            pvn_vapPressAveDailyMinSurfaceTemp = self.organicLiquidMixtureVaporMolecularWeight()['pln'] 
                 
        else:
            raise ValueError('Incorrect Assumption')
            
        pva_vaporPressureAverageDaily = pva_trueVaporPressure 
        
        # Equation 1.10.
        delta_pb = self.breatherPressureSetting - self.breatherVaccumSetting
        
        # Equation 1.9
        delta_pv = pvx_vapPressAveDailyMaxSurfaceTemp - pvn_vapPressAveDailyMinSurfaceTemp
        
        pt_4 = (deltatv_averageDailyVaporTempRange / tla_averageDailyLiquidTemp)
        pt_5 = locationTempData['atmPressure'] - pva_vaporPressureAverageDaily
        pt_6 = delta_pv - delta_pb

        ke_vaporSpaceExpansionFactor =  pt_4 + (pt_6 / pt_5)
        
        if ke_vaporSpaceExpansionFactor <= 1 and ke_vaporSpaceExpansionFactor > 0:
            _ke_vaporSpaceExpansionFactor = ke_vaporSpaceExpansionFactor
            
        else:
            raise ValueError('Assumption Invalid')
            
        return {'quantity': 'Vapor Space Expansion Factor (ke)', 
                'equation': '1-5',
                'version': '06/2020',
                'value': _ke_vaporSpaceExpansionFactor, 
                'elements': [{'delta_tv': deltatv_averageDailyVaporTempRange, 
                             'tla': tla_averageDailyLiquidTemp, 
                              'atm_pressure': locationTempData['atmPressure'], 
                              'pva': pva_vaporPressureAverageDaily,
                              'delta_pv': delta_pv, 
                              'delta_pb': delta_pb, 
                              'tlx': tlx_averageDailyMaxLiqSurfaceTemp_R, 
                              'tln': tln_averageDailyMinLiqSurfaceTemp_R, 
                              'pvx': pvx_vapPressAveDailyMaxSurfaceTemp, 
                              'pvn':pvn_vapPressAveDailyMinSurfaceTemp,
                              'product class': _productMix[0]['productClass'], 
                              'product name': _productMix[0]['productName'],
                              'status': 'None'}]}
    
    def standingLosses(self):
        '''
        Done -- Check.
        Equation 1.2.
        '''
        standingLoss = 365 *\
                        self.vaporSpaceVolume()['value'] *\
                        self.stockVaporDensity() *\
                        self.vaporSpaceExpansionFactor()['value'] *\
                        self.ventedVaporSpaceSatFactor()['value']

        return {'quantity': 'standingLoss', 
                'value': standingLoss,
                'equation': '1-2',
                'version': '06/2020',
                'elements': [{'constant': '365', 
                             'vv': self.vaporSpaceVolume()['value'] , 
                             'wv': self.stockVaporDensity(),
                             'ke': self.vaporSpaceExpansionFactor()['value'], 
                             'ks': self.ventedVaporSpaceSatFactor()['value'], 
                             'status': 'Done'}]}
  
  
productMix = [{'productClass': 'Crude Oils', 'productName': 'Crude oil (RVP 5)', 'composition': 1}]  
#productMix = [{'productClass': 'Petroleum Distillates', 'productName': 'Gasoline (RVP 8)', 'composition': 1}]  
#productMix = [{'productClass': 'Organic Liquids', 'productName': 'Benzene', 'composition': 2812/3171}, 
#              {'productClass': 'Organic Liquids', 'productName': 'Toluene', 'composition': 258/3171}, 
#              {'productClass': 'Organic Liquids', 'productName': 'Cyclohexane', 'composition': 101/3171}]
    
obj = FixedRoofTank(
              location='Long Beach, CA',
              shellShade='White', 
              shellCondition='Average',
              roofShade='White',
              roofCondition='Average',
              productMix=productMix, 
              isShellInsulated=False, 
              isRoofInsulated=False, 
              shellHeight=20,
              tankDiameter=15, 
              tankOrientation ='Vertical', 
              isAverageLiquidHeightKnown=False, 
              tankLength=20, 
              roofType='Cone',
              tankDomeRoofRadius=10, 
              roofSlope=0.0625,                  
              breatherPressureSetting=0.03, 
              breatherVaccumSetting=-0.03)
