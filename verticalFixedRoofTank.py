import math
import os 
import numpy as np
import pandas as pd
import sqlite3

class VerticalFixedRoofTank:
    
    '''
    As of now, I think the code works for single component mixtures of 
    crude oils and distillates.

    #TODO: Build out functionality for liquid mixtures, 
    and for custom data (crude and distillates)

    # TODO: Need to include handling of data that is not used, 
    but needed to instantiate the object.
    
    # TODO: Look into the use of tb for a non-insulated tank.
    '''
    
    conn = sqlite3.connect('tanks-4-09-data.db', check_same_thread=False)

    def __init__(
        self,
        shellHeight,
        tankLength,
        tankDiameter, 
        tankDomeRoofRadius, 
        isTankDomeRoofRadiusKnown,
        isTankOnVaporBalance,
        roofType, 
        tankOrientation, 
        throughput, 
        productType,
        breatherVentSetting,
        breatherPressureSetting,
        breatherVaccumSetting,
        location, 
        shellShade, 
        shellCondition, 
        roofShade, 
        roofCondition, 
        isShellInsulated,
        isRoofInsulated, 
        productClassString, 
        productNameString, 
        roofSlope, 
        isAverageLiquidHeightKnown,
        averageLiquidHeight, 
        rvp_crudeOils, 
        isCustomProduct, 
        productMW, 
        tankName, 
        assetNumber):
        
        self.hs_shellHeight = shellHeight
        self.tankDiameter = tankDiameter
        self.rr_tankDomeRoofRadius = tankDomeRoofRadius
        self.isTankDomeRoofRadiusKnown = isTankDomeRoofRadiusKnown
        self.roofType = roofType
        self.tankOrientation = tankOrientation
        self.tankLength = tankLength
        self.isTankOnVaporBalance = isTankOnVaporBalance
        self.throughput = throughput
        self.productType = productType
        self.breatherVentSetting = breatherVentSetting
        self.breatherPressureSetting = breatherPressureSetting,
        self.breatherVaccumSetting = breatherVaccumSetting,
        self.location = location
        self.shellShade = shellShade
        self.shellCondition = shellCondition
        self.roofShade = roofShade
        self.roofCondition = roofCondition
        self.isShellInsulated = isShellInsulated
        self.isRoofInsulated = isRoofInsulated
        self.productClassString = productClassString
        self.productNameString = productNameString
        self.roofSlope = roofSlope
        self.isAverageLiquidHeightKnown = isAverageLiquidHeightKnown
        self.averageLiquidHeight = averageLiquidHeight
        self.rvp_crudeOils = rvp_crudeOils
        self.isCustomProduct = isCustomProduct
        self.productMW = productMW
        self.tankName = tankName
        self.assetNumber = assetNumber
        
        # Get data from db for the location.
        locationString1 = self.location
        locationString2 = "SELECT * FROM revised_met_data WHERE location == "
        locationString3 = "'" + str(locationString1) + "'"
        queryStringTax = locationString2 + locationString3 
        self.sqlQuery_ = pd.read_sql(queryStringTax, self.conn)
        
        # Parse out TAX, TAN, wind speed, solar data, and atm pressure from the sql query.
        # Assume that the TAX and TAN values are the annual values. 
        self.tax_maxAmbientTemp_R = self.sqlQuery_[self.sqlQuery_['symbol'] == 'TAX'].values.tolist()[0][-1] + 460.67
        self.tan_minAmbientTemp_R = self.sqlQuery_[self.sqlQuery_['symbol'] == 'TAN'].values.tolist()[0][-1] + 460.67
        
        if self.productClassString in ['Crude Oils']:
            if self.productNameString in ['Crude oil (RVP 5)']:
        
                _productClassString = self.productClassString
                _productNameString = self.productNameString
                _queryStringProduct1 = "SELECT * FROM chemical_data " 
                _queryStringProduct2 = "WHERE CATEGORY == " + "'"  
                _queryStringProduct3 = str(productClassString) + "'"
                _queryStringProduct = _queryStringProduct1 + _queryStringProduct2 + _queryStringProduct3
                _sqlQueryProduct = pd.read_sql(_queryStringProduct, self.conn)
                _sqlQueryProductData = _sqlQueryProduct[_sqlQueryProduct['NAME'] == _productNameString]
                self.sqlQueryProductRVP = _sqlQueryProductData['REID'].values.tolist()[0]
                self.sqlQueryProductMW = _sqlQueryProductData['VP_MOLWT'].values.tolist()[0]
                self.sqlQueryProductNAME = _sqlQueryProductData['NAME'].values.tolist()[0]
            
            else:
                raise ValueError('Assumption Invalid')
                
        elif self.productClassString in ['Petroleum Distillates']:
            if self.productNameString in ['Gasoline (RVP 10)', 'Gasoline (RVP 13)', 'Gasoline (RVP 7)',
                                          'Gasoline (RVP 6)', 'Gasoline (RVP 8)', 'Gasoline (RVP 9)',
                                          'Gasoline (RVP 11)', 'Gasoline (RVP 12)', 'Gasoline (RVP 7.8)',
                                          'Gasoline (RVP 8.3)', 'Gasoline (RVP 11.5)','Gasoline (RVP 13.5)', 
                                          'Gasoline (RVP 15.0)']:
                
                _productClassString = self.productClassString
                _productNameString = self.productNameString
                _queryStringProduct1 = "SELECT * FROM chemical_data " 
                _queryStringProduct2 = "WHERE CATEGORY == " + "'"  
                _queryStringProduct3 = str(productClassString) + "'"
                _queryStringProduct = _queryStringProduct1 + _queryStringProduct2 + _queryStringProduct3
                _sqlQueryProduct = pd.read_sql(_queryStringProduct, self.conn)
                _sqlQueryProductData = _sqlQueryProduct[_sqlQueryProduct['NAME'] == _productNameString]
                self.sqlQueryProductRVP = _sqlQueryProductData['REID'].values.tolist()[0]
                self.sqlQueryProductRVP = _sqlQueryProductData['REID'].values.tolist()[0]
                self.sqlQueryProductMW = _sqlQueryProductData['VP_MOLWT'].values.tolist()[0]
                self.sqlQueryProductNAME = _sqlQueryProductData['NAME'].values.tolist()[0]
            
            else:
                raise ValueError('Assumption Invalid')
        
        elif self.productClassString in ['Organic Liquids']:
            pass
        
        else:
            raise ValueError('Assumption Invalid') 

        
    def averageAmbientTempData(self):

        ''' 
        Partially done -- Look into tb for fully insulated 
                          or partially insulated tanks.
        Equation 1-11.
        '''
        
        data_list_wind = self.sqlQuery_[self.sqlQuery_['symbol'] == 'V'].values.tolist()[0][-1]
        i_solarInsulation = self.sqlQuery_[self.sqlQuery_['symbol'] == 'I'].values.tolist()[0][-1]
        atmPressure = self.sqlQuery_[self.sqlQuery_['symbol'] == 'PA'].values.tolist()[0][-1]
        
        # Equation 1.11. Assumed to be the difference between the annual value.
        delta_ta_R = (self.tax_maxAmbientTemp_R - self.tan_minAmbientTemp_R)
        
        # Equation 1.30, value is in Rankine. Assumed to be the sum 
        #    of the annual values divided by 2
        taa_averageDailyTemp = (self.tan_minAmbientTemp_R + self.tax_maxAmbientTemp_R) * (1/2.)

        # Get data from for the shell/roof condition and color.
        queryStringShell1 = " FROM solar_abs_data WHERE Surface_Color == " + "'" + str(self.shellShade) + "'"
        queryStringShell = "SELECT " + str(self.shellCondition) + queryStringShell1
        sqlQueryqueryStringShell = pd.read_sql(queryStringShell, self.conn)
        alphaS_shell = sqlQueryqueryStringShell.iloc[0][self.shellCondition]

        queryStringRoof1 = " FROM solar_abs_data WHERE Surface_Color == " + "'" + str(self.roofShade) + "'"
        queryStringRoof = "SELECT " + str(self.roofCondition) + queryStringRoof1
        sqlQueryqueryStringRoof = pd.read_sql(queryStringRoof, self.conn)
        alphaR_roof = sqlQueryqueryStringRoof.iloc[0][self.roofCondition]
        
        # Equation 1.31. Assumes that the tank is not insulated. 
        # Not sure about tb for fully insulated or partially insulated tanks.
        tb_liquidBulkTemp = taa_averageDailyTemp + (0.003*alphaS_shell*i_solarInsulation)

        return {'city': self.location,
                'shell_shade': self.shellShade,
                'shell_condition': self.shellCondition,
                'roof_shade': self.roofShade,
                'roof_condition': self.roofCondition,
                'tax_maxAmbientTemp_R': self.tax_maxAmbientTemp_R, 
                'tan_minAmbientTemp_R': self.tan_minAmbientTemp_R, 
                'wind_speed': data_list_wind,
                'delta_ta_R':delta_ta_R, 
                'alphaS_shell': alphaS_shell, 
                'alphaR_roof': alphaR_roof, 
                'i_solarInsulation': i_solarInsulation,
                'taa_averageDailyTemp': taa_averageDailyTemp, 
                'atmPressure': atmPressure, 
                'tb_liquidBulkTemp': tb_liquidBulkTemp, 
                'version': '06/2020, Table: 7.1-7', 
                'status': 'Partially Done'}
    
    def vaporSpaceVolume(self):
        '''
        Done - check.
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

    # ========================================================================= #
    
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
        pva_vaporPressureAverageDaily = self.vaporSpaceExpansionFactor()['elements'][0]['pva']   
        
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
    
    # ========================================================================= #
    
    def vaporSpaceOutage(self):
        '''
        Done - Check
        Equation 1-16. 
        Uses Notes 1, 2 and Equations 1-17 - 1-20.
        '''

        # Assumptions made within the function.
        rs_tankShellRadius = self.tankDiameter / 2.
        
        if rs_tankShellRadius >=0.8*self.tankDiameter or\
            rs_tankShellRadius <= 1.2*self.tankDiameter:
            pass
        
        # Assumed to be half of the shell height.
        if self.isAverageLiquidHeightKnown == True:
            hl_liquidLevel = self.averageLiquidHeight
            
        elif self.isAverageLiquidHeightKnown == False:
            hl_liquidLevel = self.hs_shellHeight / 2.
            
        else:
            raise ValueError('Invalid Assumption')

        # hro_roof_outage
        if self.roofType == 'Dome':

            if self.isTankDomeRoofRadiusKnown == True:
                
                valueSqRoot = (self.rr_tankDomeRoofRadius**2) -\
                              (rs_tankShellRadius**2)
                
                if valueSqRoot >= 0:
                    #Equation 1.20.
                    hr_roofHeight   =   self.rr_tankDomeRoofRadius -\
                                        (valueSqRoot)**(0.5)
                
                else:
                    raise ValueError('Taking the negative of a square root')
                
                # Equation 1.19
                hro_roofOutage  =   hr_roofHeight *\
                                    ( (1 / 2.) + ( (1 / 6.) *\
                                    (hr_roofHeight / rs_tankShellRadius)**2) )
            
            elif self.isTankDomeRoofRadiusKnown == False:
                
                #Why does this matter if the value is unknown?
                hr_roofHeight = 0.268 * rs_tankShellRadius 
                hro_roofOutage = 0.137 * rs_tankShellRadius
            
            else:
                raise ValueError('Assumption Invalid')

        elif self.roofType == 'Cone':
            # Note 1 to determine hro_roof_outage for a cone roof.
            # Equation 1.17, and 1.18.
            hro_roofOutage = ( 1 / 3. ) * self.roofSlope * rs_tankShellRadius
        
        else:
            raise ValueError('Incorrect Roof Type.')
        
        # hvo_vaporSpaceOutage
        if self.tankOrientation == 'Vertical':
            #Equation 1.16.
            hvo_vaporSpaceOutage =    self.hs_shellHeight - \
                                        hl_liquidLevel + hro_roofOutage
            
            result = {'value': hvo_vaporSpaceOutage, 
                      'quantity': 'Vapor Space Outage (hvo)',
                      'equation': '1-16',
                      'version': '06/2020', 'elements': [{'hs': self.hs_shellHeight, 
                                                          'hl': hl_liquidLevel, 
                                                          'hro': hro_roofOutage, 
                                                          'hr': self.roofSlope * rs_tankShellRadius, 
                                                          'roof slope': self.roofSlope, 
                                                          'roof height': self.roofSlope * rs_tankShellRadius,
                                                          'shell radius': rs_tankShellRadius,
                                                          'shell height': self.hs_shellHeight,
                                                          'liquid height': hl_liquidLevel,
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
    
    # ========================================================================= #
    
    def stockVaporDensity(self):
        '''
        Done -- Check.
        Equation 1-22.
        '''
        
        locationTempData = self.averageAmbientTempData()
        hs_d = self.hs_shellHeight / self.tankDiameter
        
        # Fully Insulated
        if self.isShellInsulated == True and self.isRoofInsulated == True:
            
            tv_averageVaporTemp = locationTempData['tb_liquidBulkTemp']
        
        #Partially Insulated
        elif self.isShellInsulated == True and self.isRoofInsulated == False:
            
            pt_1 = 0.6*locationTempData['taa_averageDailyTemp']
            pt_2 = 0.4*locationTempData['tb_liquidBulkTemp'] 
            pt_3 = 0.01*locationTempData['alphaR_roof']*locationTempData['i_solarInsulation']
            
            pt_4 = ((2.2 * hs_d) + 1.1)*locationTempData['taa_averageDailyTemp']
            pt_5 = 0.8*locationTempData['tb_liquidBulkTemp'] 
            pt_6 = 0.021*locationTempData['alphaR_roof']*locationTempData['i_solarInsulation']
            pt_7 = 0.013*hs_d*locationTempData['alphaS_shell']*locationTempData['i_solarInsulation']  
            
            tv_averageVaporTemp = (pt_4 + pt_5 + pt_6 + pt_7) / (pt_4)
        
        #Uninsulated
        elif self.isShellInsulated == False and self.isRoofInsulated == False:
            
            tv_averageVaporTemp = locationTempData['tb_liquidBulkTemp']
            
        if self.isCustomProduct == True:
            
            mv_vaporMW = self.productMW
            name_material = self.productNameString
        
        elif self.isCustomProduct == False:
                
            queryStringProduct1 = "'" + str(self.productClassString) + "'"
            queryStringProduct = "SELECT * FROM chemical_data WHERE CATEGORY == " + queryStringProduct1

            sqlQueryProduct = pd.read_sql(queryStringProduct, self.conn)
            mv_vaporMW = self.sqlQueryProductMW 
            name_material = self.sqlQueryProductNAME 
        
        else:
            raise ValueError('Error')

        pva_vaporPressureAverageDaily = self.vaporSpaceExpansionFactor()['elements'][0]['pva']
        
        wv_vaporDensity = (pva_vaporPressureAverageDaily * mv_vaporMW) / (10.731 * tv_averageVaporTemp)
        
        return {'quantity': 'Stock Vapor Density (wv)', 
                'value': wv_vaporDensity, 
                'equation': '1-22',
                'version': '06/2020',
                'elements': [{'pva': pva_vaporPressureAverageDaily, 
                              'name': name_material,
                              'mv': mv_vaporMW, 
                              'constant': 10.731, 
                              'tv': tv_averageVaporTemp, 
                              'status': 'Done'}]}
    
    # ========================================================================= #
    
    def vaporSpaceExpansionFactor(self):
        '''
        Partially done -- Look into the value for tb 
        and the tank insulation, and the partial pressures. 
        Equation 1-5.
        
        Most of the values match the OK printout.
        
        Also, look into `delta_pv`, and need to fold in 
            logifc to pull in rvp based in user input.
    
        '''
        locationTempData = self.averageAmbientTempData()
                
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
            hs_d = self.hs_shellHeight / self.tankDiameter
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
                
        # Only for Crude Oils, and Selected Petroleum Stocks
        if self.productClassString in ['Crude Oils', 'Petroleum Distillates']:
            
            if self.isCustomProduct == False:
                rvp_crude = self.sqlQueryProductRVP
                
            elif self.isCustomProduct == True:
                rvp_crude = self.rvp_crudeOils
                
            else:
                raise ValueError('Invalid Assumption')
                
            const_a = 12.82 - 0.9672*math.log(rvp_crude)
            const_b = 7261 - 1216*math.log(rvp_crude)
            
            constDict =  {'const_a': const_a, 'const_b': const_b}
            val_pva = constDict['const_a'] - ( (constDict['const_b']) / tla_averageDailyLiquidTemp )
            val_pvx = constDict['const_a'] - ( (constDict['const_b']) / tlx_averageDailyMaxLiqSurfaceTemp_R )
            val_pvn = constDict['const_a'] - ( (constDict['const_b']) / tln_averageDailyMinLiqSurfaceTemp_R )
            
            pva_trueVaporPressure = math.exp(val_pva)
            pvx_vapPressAveDailyMaxSurfaceTemp = math.exp(val_pvx)
            pvn_vapPressAveDailyMinSurfaceTemp = math.exp(val_pvn)
        
        # TODO: Need to resolve this function for organic liquids.
        elif self.productClassString in ['Organic Liquids']:
            # Use Equation 1.26. Need to get data from a higher level due to scoping.
            pass
                 
        
        else:
            pass
            raise ValueError('Incorrect Assumption')
            
        pva_vaporPressureAverageDaily = pva_trueVaporPressure 
        
        # Equation 1.10.
        delta_pb = self.breatherPressureSetting[0] - self.breatherVaccumSetting[0]
        
        # Equation 1.9
        delta_pv = pvx_vapPressAveDailyMaxSurfaceTemp - pvn_vapPressAveDailyMinSurfaceTemp
        
        pt_4 = (deltatv_averageDailyVaporTempRange / tla_averageDailyLiquidTemp)
        pt_5 = locationTempData['atmPressure'] - pva_vaporPressureAverageDaily
        pt_6 = delta_pv - delta_pb

        ke_vaporSpaceExpansionFactor =  pt_4 + (pt_6/pt_5)

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
                              'const_a': const_a,
                              'const_b': const_b,
                              'delta_pv': delta_pv, 
                              'delta_pb': delta_pb, 
                              'tlx': tlx_averageDailyMaxLiqSurfaceTemp_R, 
                              'tln': tln_averageDailyMinLiqSurfaceTemp_R, 
                              'pvx': pvx_vapPressAveDailyMaxSurfaceTemp, 
                              'pvn':pvn_vapPressAveDailyMinSurfaceTemp,
                              'rvp': rvp_crude,
                              'product class': self.productClassString, 
                              'product name': self.productNameString,
                              'status': 'Partially Done. Need to look into tb \
                                         for uninsulated tanks, and fold in \
                                         the vapor pressures for org. liquids'}]}
    
    # ========================================================================= #
    
    def totalLosses(self):
        '''
        Done -- check.
        Equation 1.1.
        '''
        
        totalLoss = self.workingLosses()['value'] + self.standingLosses()['value']
        
        return {'quantity': 'totalLosses', 
                'equation': '1-1',
                'version': '06/2020',
                'value': totalLoss, 
                'elements': [{'lw': self.workingLosses()['value'],
                              'ls': self.standingLosses()['value'], 
                              'status': 'Done'}]}
    # ========================================================================= #
    
    def standingLosses(self):
        '''
        Done -- Check.
        Equation 1.2.
        '''
        standingLoss = 365 *\
                        self.vaporSpaceVolume()['value']*\
                        self.stockVaporDensity()['value']*\
                        self.vaporSpaceExpansionFactor()['value']*\
                        self.ventedVaporSpaceSatFactor()['value']

        return {'quantity': 'standingLoss', 
                'value': standingLoss,
                'equation': '1-2',
                'version': '06/2020',
                'elements': [{'constant': '365', 
                             'vv': self.vaporSpaceVolume()['value'] , 
                             'wv': self.stockVaporDensity()['value'],
                             'ke': self.vaporSpaceExpansionFactor()['value'], 
                             'ks': self.ventedVaporSpaceSatFactor()['value'], 
                             'status': 'Done'}]}
    # ========================================================================= #
    
    def workingLosses(self):

        '''
        Done -- Check.
        Equation 1-35.
        '''
        pa_atmosphericPressure = self.averageAmbientTempData()['atmPressure']
        # ========================================================================= #
        
        # hqi_annualSumIncreases
        if self.tankOrientation == 'Horizontal':
            effectiveDiameter = ( (self.tankLength * self.tankDiameter) / \
                                 ( (1. / 4.) * (math.pi) ) ) ** (0.5)
            _tankDiameter = effectiveDiameter
        
        elif self.tankOrientation == 'Vertical':
            _tankDiameter = self.tankDiameter
        
        else:
            raise ValueError('Assumption Invalid')
        
        hqi_annualSumIncreases = (5.614 * self.throughput/42) /\
                                    (((math.pi) / 4.) * (_tankDiameter)**(2))
         # ========================================================================= #       
        
        # hlx_maximumLiquidHeight, and hln_minimumLiquidHeight
        # Assume that in both cases, the min. and max liquid level is not known. 
        if self.tankOrientation == 'Horizontal':
            # In this case, `self.tankDiameter` is the vertical cross section of the horizontal tank.
            hlx_maximumLiquidHeight = (math.pi) * (1. / 4.) * self.tankDiameter
            hln_minimumLiquidHeight = 0

        elif self.tankOrientation == 'Vertical':
            
            hlx_maximumLiquidHeight = self.hs_shellHeight - 1
            hln_minimumLiquidHeight = 1

        else:
            raise ValueError('Assumption Invalid')
        
        # ========================================================================= #      
         
        # n_numberTurnovers
        n_numberTurnovers = hqi_annualSumIncreases /\
                            (hlx_maximumLiquidHeight - hln_minimumLiquidHeight)
        
        # ========================================================================= #
        
        # Assume no splash loading.
        # TODO: Address tanks where flashing occurs.
        # kn_workingLossTurnOverFactor
        if self.isTankOnVaporBalance == True:
            kn_workingLossTurnOverFactor = 1
        
        elif  self.isTankOnVaporBalance == False:
            if n_numberTurnovers <= 36:
                kn_workingLossTurnOverFactor = 1

            elif n_numberTurnovers > 36:
                kn_workingLossTurnOverFactor = (180 + n_numberTurnovers) /\
                                                (6 * n_numberTurnovers)
                
            else:
                raise ValueError('Assumption Invalid')
        else:
            raise ValueError('Assumption Invalid')
            
        # ========================================================================= #
        
        # kp_workingLossProductFactor
        # Assume no splash loading into tank.
        if self.productType == 'Organic Liquids':
            kp_workingLossProductFactor = 1

        elif self.productType == 'Crude Oils':
            kp_workingLossProductFactor = 0.75

        else:
            raise ValueError('Assumption Invalid')
             
        # ========================================================================= #    
            
        # vq_netWorkingLossThroughput
        # Use Equation 1-39, check for assumption.
        if self.tankOrientation == 'Horizontal':
                
            effectiveDiameter = ( (self.tankLength * self.tankDiameter) / \
                                 ( (1. / 4.) * (math.pi) ) ) ** (0.5)
            _tankDiameter = effectiveDiameter
            
            vq_netWorkingLossThroughput = hqi_annualSumIncreases * (math.pi/4.)*(_tankDiameter)**2
        
        elif self.tankOrientation == 'Vertical':
            
            _tankDiameter = self.tankDiameter
            vq_netWorkingLossThroughput = hqi_annualSumIncreases * (math.pi/4.)*((_tankDiameter)**2)

        else:
            raise ValueError('Assumption Invalid')
    
        # ========================================================================= #
   
        pva_vaporPressureAverageDaily = self.vaporSpaceExpansionFactor()['elements'][0]['pva']

        # TODO: Look at Custom Case more carefully. 
        if self.breatherVentSetting == 'Custom':
            
            pbp_breatherVentPressureSetting = self.breatherPressureSetting
            
            _condition =    kn_workingLossTurnOverFactor * \
                            ((pbp_breatherVentPressureSetting + pa_atmosphericPressure) / \
                            (pi_pressureVaporSpace + pa_atmosphericPressure))
            
            if _condition > 1:
                
                kb_ventSettingCorrection = (((pi_pressureVaporSpace + pa_atmosphericPressure) / \
                             kn_workingLossTurnOverFactor) - pva_vaporPressureAverageDaily) / \
                            (pbp_breatherVentPressureSetting + pa_atmosphericPressure - \
                             pva_vaporPressureAverageDaily)
            
            elif _condition <= 1:
                kb_ventSettingCorrection=  kb_ventSettingCorrection
            
            else:
                raise ValueError('Assumption')
            
        # Equation 1-41
        elif self.breatherVentSetting == 'Default':
            kb_ventSettingCorrection = 1 #self.breatherPressureSetting[0] 
        
        else:
            raise ValueError('Assumption Invalid')
        
        workingLoss =   vq_netWorkingLossThroughput *\
                        kn_workingLossTurnOverFactor *\
                        kp_workingLossProductFactor *\
                        self.stockVaporDensity()['value'] *\
                        kb_ventSettingCorrection
        
        return {'quantity': 'workingLoss', 
                'value': round(workingLoss, 4),
                'equation': '1-2',
                'version': '06/2020',
                'elements': [{'tankOrientation': self.tankOrientation,
                              'vent setting': self.breatherVentSetting,
                              'pva': pva_vaporPressureAverageDaily,
                             'vq': vq_netWorkingLossThroughput, 
                             'kn': kn_workingLossTurnOverFactor, 
                             'kp': kp_workingLossProductFactor,
                             'wv': self.stockVaporDensity()['value'], 
                             'kb': kb_ventSettingCorrection, 
                              'hlx': hlx_maximumLiquidHeight, 
                              'hln': hln_minimumLiquidHeight, 
                              'q_gal': self.throughput, 
                              'q_bbl': self.throughput / 42,
                              'diameter': _tankDiameter,
                              'num_turnovers': n_numberTurnovers, 
                              'sum_hqi': hqi_annualSumIncreases,
                              'kbv': kb_ventSettingCorrection,
                              'kbp': self.breatherVaccumSetting[0],
                              'on_vapor_bal':  self.isTankOnVaporBalance,
                              'breatherPressureSetting': self.breatherPressureSetting, 
                              'breatherVaccumSetting': self.breatherVaccumSetting,
                              'roof insulated': self.isRoofInsulated,
                              'shell insulated': self.isShellInsulated,
                             'status': 'Done'}]}
    
    def tankInputs(self):
        
        dataDict = {
         'shellHeight': self.hs_shellHeight,
         'tankLength': self.tankLength,
         'tankDiameter': self.tankDiameter, 
         'tankDomeRoofRadius': self.rr_tankDomeRoofRadius, 
         'isTankDomeRoofRadiusKnown': self.isTankDomeRoofRadiusKnown,
         'isTankOnVaporBalance': self.isTankOnVaporBalance,
         'roofType': self.roofType, 
         'tankOrientation': self.tankOrientation, 
         'throughput': self.throughput, 
         'productType': self.productType,
         'breatherVentSetting': self.breatherVentSetting,
         'breatherPressureSetting': self.breatherPressureSetting,
         'breatherVaccumSetting': self.breatherVaccumSetting,
         'location': self.location, 
         'shellShade': self.shellShade, 
         'shellCondition': self.shellCondition, 
         'roofShade': self.roofShade, 
         'roofCondition': self.roofCondition, 
         'isShellInsulated': self.isShellInsulated,
         'isRoofInsulated': self.isRoofInsulated, 
         'productClassString': self.productClassString, 
         'productNameString': self.productNameString, 
         'roofSlope': self.roofSlope, 
         'isAverageLiquidHeightKnown': self.isAverageLiquidHeightKnown,
         'averageLiquidHeight': self.averageLiquidHeight, 
         'rvp_crudeOils': self.rvp_crudeOils, 
         'tankName': self.tankName,
         'assetNumber': self.assetNumber}

        return dataDict
