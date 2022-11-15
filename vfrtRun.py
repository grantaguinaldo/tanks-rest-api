import json
import math
from sqlite3 import connect
import sqlite3
import pandas as pd
from verticalFixedRoofTank import VerticalFixedRoofTank

def returnData1(shellHeight=None, 
            tankDiameter=None, 
            tankDomeRoofRadius=None, 
            breatherPressureSetting=None,
            breatherVaccumSetting=None,
            tankLength=None,
            throughput=None, 
            roofSlope=None, 
            averageLiquidHeight=None, 
            rvp_crudeOils=None,
            productMW=None,
            isTankDomeRoofRadiusKnown=None, 
            isTankOnVaporBalance=None, 
            isRoofInsulated=None,
            isShellInsulated=None,   
            isAverageLiquidHeightKnown=None,
            isCustomProduct=None,
            roofType=None,
            tankOrientation=None, 
            productType=None, 
            breatherVentSetting=None,
            location=None,
            shellShade=None,
            shellCondition=None, 
            roofShade=None,
            roofCondition=None, 
            productClassString=None,
            productNameString=None,         
            tankName=None,
            assetNumber=None):
        
    obj =   VerticalFixedRoofTank(
            shellHeight=shellHeight, 
            tankDiameter=tankDiameter, 
            tankDomeRoofRadius=tankDomeRoofRadius, 
            breatherPressureSetting=breatherPressureSetting,
            breatherVaccumSetting=breatherVaccumSetting,
            tankLength=tankLength,
            throughput=throughput, 
            roofSlope=roofSlope, 
            averageLiquidHeight=averageLiquidHeight, 
            rvp_crudeOils=rvp_crudeOils,
            productMW=productMW,
            isTankDomeRoofRadiusKnown=isTankDomeRoofRadiusKnown, 
            isTankOnVaporBalance=isTankOnVaporBalance, 
            isRoofInsulated=isRoofInsulated,
            isShellInsulated=isShellInsulated,   
            isAverageLiquidHeightKnown=isAverageLiquidHeightKnown,
            isCustomProduct=isCustomProduct,
            roofType=roofType,
            tankOrientation=tankOrientation, 
            productType=productType, 
            breatherVentSetting=breatherVentSetting,
            location=location,
            shellShade=shellShade,
            shellCondition=shellCondition, 
            roofShade=roofShade,
            roofCondition=roofCondition, 
            productClassString=productClassString,
            productNameString=productNameString,         
            tankName=tankName,
            assetNumber=assetNumber)

    a = obj.averageAmbientTempData()
    b = obj.vaporSpaceVolume()
    c = obj.ventedVaporSpaceSatFactor()
    d = obj.vaporSpaceOutage()
    e = obj.stockVaporDensity()
    f = obj.vaporSpaceExpansionFactor()
    g = obj.totalLosses()
    h = obj.standingLosses()
    i = obj.workingLosses()
    j = obj.tankInputs()

    data = [ {'index': 0, 'method': 'averageAmbientTempData', 'data': a}, 
             {'index': 1, 'method': 'vaporSpaceVolume', 'data': b}, 
             {'index': 2, 'method': 'ventedVaporSpaceSatFactor', 'data': c}, 
             {'index': 3, 'method': 'vaporSpaceOutage', 'data': d},
             {'index': 4, 'method': 'stockVaporDensity', 'data': e}, 
             {'index': 5, 'method': 'vaporSpaceExpansionFactor', 'data': f}, 
             {'index': 6, 'method': 'totalLosses', 'data': g}, 
             {'index': 7, 'method': 'standingLosses', 'data': h},
             {'index': 8, 'method': 'workingLosses', 'data': i}, 
             {'index': 9, 'method': 'tankInputs', 'data': j}]
    
    #tank dimensions
    hs = data[9]['data']['shellHeight']
    d = data[8]['data']['elements'][0]['diameter']
    hlx = data[8]['data']['elements'][0]['hlx']
    hln = data[8]['data']['elements'][0]['hln']
    hl = data[3]['data']['elements'][0]['hl'] #ok
    n = data[8]['data']['elements'][0]['num_turnovers']
    throughput_gal = data[8]['data']['elements'][0]['q_gal']
    q_bbl = data[8]['data']['elements'][0]['q_bbl']
    on_vapor_bal = data[8]['data']['elements'][0]['on_vapor_bal']
    
    #paint characteristics
    shell_color = data[0]['data']['shell_shade']
    shell_cond = data[0]['data']['shell_condition']
    roof_color = data[0]['data']['roof_shade']
    roof_cond = data[0]['data']['roof_condition']
    
    #roof characteristics
    roof_type = data[9]['data']['roofType']
    hr = data[3]['data']['elements'][0]['hr']
    slope = data[9]['data']['roofSlope']
    
    #Breather Vent
    pbv = data[8]['data']['elements'][0]['breatherVaccumSetting']
    pbp = data[8]['data']['elements'][0]['breatherPressureSetting']
    
    #Insulation
    is_shell_insulated = data[8]['data']['elements'][0]['shell insulated']
    is_roof_insulated = data[8]['data']['elements'][0]['roof insulated']
    
    #MET Data
    city = data[0]['data']['city']
    taa = data[0]['data']['taa_averageDailyTemp']
    tax = data[0]['data']['tax_maxAmbientTemp_R'] 
    tan = data[0]['data']['tan_minAmbientTemp_R']
    del_ta = data[0]['data']['delta_ta_R']
    v =  data[0]['data']['wind_speed'] 
    i = data[0]['data']['i_solarInsulation'] 
    pa = data[0]['data']['atmPressure'] 
    
    #Liquid Data
    liq_category = data[5]['data']['elements'][0]['product class']
    liq_name = data[5]['data']['elements'][0]['product name']
    tb = data[0]['data']['tb_liquidBulkTemp']
    tla = data[5]['data']['elements'][0]['tla']
    tln = data[5]['data']['elements'][0]['tln']
    tlx = data[5]['data']['elements'][0]['tlx']
    pva = data[5]['data']['elements'][0]['pva']
    pvn = data[5]['data']['elements'][0]['pvn']
    pvx = data[5]['data']['elements'][0]['pvx']
    mv = data[4]['data']['elements'][0]['mv']
    rvp = data[5]['data']['elements'][0]['rvp']
    const_a = data[5]['data']['elements'][0]['const_a']
    const_b = data[5]['data']['elements'][0]['const_b']
    
    #Standing Losses
    ls = data[7]['data']['value']
    vv = data[7]['data']['elements'][0]['vv']
    wv= data[7]['data']['elements'][0]['wv']
    ke= data[7]['data']['elements'][0]['ke']
    ks= data[7]['data']['elements'][0]['ks']
    
    #Working Losses
    lw = data[8]['data']['value']
    vq = data[8]['data']['elements'][0]['vq']   
    kn = data[8]['data']['elements'][0]['kn']   
    kp = data[8]['data']['elements'][0]['kp']
    wv = data[8]['data']['elements'][0]['wv']   
    kb = data[8]['data']['elements'][0]['kb']
    
    #Total Losses
    lw = data[6]['data']['elements'][0]['lw'] 
    ls = data[6]['data']['elements'][0]['ls'] 
    lt = data[6]['data']['value']
    
    #Vapor Space Volume
    vv = data[7]['data']['elements'][0]['vv']
    diameter = data[8]['data']['elements'][0]['diameter']  
    hvo = data[1]['data']['elements'][0]['hvo']
    
    #Vapor Space Outage
    hvo = data[1]['data']['elements'][0]['hvo']
    hs = data[9]['data']['shellHeight']
    hl = data[3]['data']['elements'][0]['hl'] #ok
    hro = data[3]['data']['elements'][0]['hro'] #ok
    
    #Roof Outage
    hro = data[3]['data']['elements'][0]['hro'] #ok
    hr = data[3]['data']['elements'][0]['hr']
    rs = data[9]['data']['tankDiameter'] / 2.
    sr = data[9]['data']['roofSlope']
    
    #Vapor Density
    wv = data[8]['data']['elements'][0]['wv']
    mv = data[4]['data']['elements'][0]['mv']
    pva = data[5]['data']['elements'][0]['pva']
    r = data[4]['data']['elements'][0]['constant']
    tv = data[4]['data']['elements'][0]['tv']
    alpha_s = data[0]['data']['alphaS_shell']
    alpha_r = data[0]['data']['alphaR_roof']
    i_value = data[0]['data']['i_solarInsulation']
    
    #Vapor Space Expansion Factor
    ke= data[7]['data']['elements'][0]['ke']
    del_tv = data[5]['data']['elements'][0]['delta_tv']
    delta_pv = data[5]['data']['elements'][0]['delta_pv']
    delta_pb = data[5]['data']['elements'][0]['delta_pb']
    pa = data[5]['data']['elements'][0]['atm_pressure']
    pva = data[5]['data']['elements'][0]['pva']
    tla = data[5]['data']['elements'][0]['tla']

    #Vented Vapor Space Saturation Factor
    ks = data[2]['data']['value']
    pva = data[2]['data']['elements'][0]['pva']
    hvo = data[2]['data']['elements'][0]['hvo']
    
    return_data = {

        'tank_dimensions': [{

            'value': None,
            'parameter_estimates': [{'parameter_name': 'Shell Height, ft', 'parameter_symbol': 'hs', 'parameter_value': hs},
                                    {'parameter_name': 'Shell Diameter, ft', 'parameter_symbol': 'd', 'parameter_value': d},
                                    {'parameter_name': 'Max Liquid Height, ft', 'parameter_symbol': 'hlx', 'parameter_value': hlx},
                                    {'parameter_name': 'Min Liquid Height, ft', 'parameter_symbol': 'hln', 'parameter_value': hln},
                                    {'parameter_name': 'Average Liquid Height, ft', 'parameter_symbol': 'hl', 'parameter_value': hl},
                                    {'parameter_name': 'Number of Turnovers per Year, unitless', 'parameter_symbol': 'n', 'parameter_value': n},
                                    {'parameter_name': 'Net Annual Throughput, gal', 'parameter_symbol': 'q_gal', 'parameter_value': throughput_gal},
                                    {'parameter_name': 'Net Annual Throughput, bbl', 'parameter_symbol': 'q_bbl', 'parameter_value': q_bbl},
                                    {'parameter_name': 'Vapor Balanced, bool', 'parameter_symbol': 'on_vapor_bal','parameter_value': on_vapor_bal}]
            }],

        'paint_characteristics': [{

            'value': None,
            'parameter_estimates': [{'parameter_name': 'Shell Color', 'parameter_symbol': 'shell_color', 'parameter_value': shell_color},
                                    {'parameter_name': 'Shell Condition', 'parameter_symbol': 'shell_cond', 'parameter_value': shell_cond},
                                    {'parameter_name': 'Roof Color', 'parameter_symbol': 'roof_color', 'parameter_value': roof_color},
                                    {'parameter_name': 'Roof Condition', 'parameter_symbol': 'roof_cond', 'parameter_value': roof_cond}]
        }],

        'roof_characteristics': [{

            'value': None,
            'parameter_estimates': [{'parameter_name': 'roof type', 'parameter_symbol': 'roof_type', 'parameter_value': roof_type},
                                    {'parameter_name': 'roof height, ft', 'parameter_symbol': 'hr', 'parameter_value': hr},
                                    {'parameter_name': 'roof slope, ft/ft', 'parameter_symbol': 'slope', 'parameter_value': slope}]
        }],


        'breather_vent_settings': [{

            'value': None,
            'parameter_estimates': [{'parameter_name': 'Breather Vent Vacuum Setting', 'parameter_symbol': 'pbv, psig', 'parameter_value': pbv[0]},
                                    {'parameter_name': 'Breather Vent Pressure Setting', 'parameter_symbol': 'pbp, psig', 'parameter_value': pbp[0]}]
        }],


        'insulation_characteristics': [{

            'value': None,
            'parameter_estimates': [{'parameter_name': 'Shell Insulation', 'parameter_symbol': 'is_shell_insulated, bool', 'parameter_value': is_shell_insulated},
                                    {'parameter_name': 'Roof Insulation', 'parameter_symbol': 'is_roof_insulated, bool', 'parameter_value': is_roof_insulated}]
        }],

        'meterological_data': [{

            'value': None,
            'parameter_estimates': [{'parameter_name': 'City', 'parameter_symbol': 'city', 'parameter_value': city},
                                    {'parameter_name': 'Average Daily Ambient Temp', 'parameter_symbol': 'taa, R', 'parameter_value': taa},
                                    {'parameter_name': 'Average Daily Max Ambient Temp', 'parameter_symbol': 'tax, R', 'parameter_value': tax},
                                    {'parameter_name': 'Average Daily Min Ambient Temp', 'parameter_symbol': 'tan, R', 'parameter_value': tan},
                                    {'parameter_name': 'Average Daily Ambient Temp Range', 'parameter_symbol': 'del_ta, R', 'parameter_value': del_ta},
                                    {'parameter_name': 'Average Wind Speed', 'parameter_symbol': 'v, miles/hr', 'parameter_value': v},
                                    {'parameter_name': 'Average Daily Total Insolation Factor', 'parameter_symbol': 'i, btu/ft^2-day', 'parameter_value': i},
                                    {'parameter_name': 'Atm. Pressure', 'parameter_symbol': 'pa, psi', 'parameter_value': pa}]
        }],

        'liquid_data': [{

            'value': None,
            'parameter_estimates': [{'parameter_name': 'Liquid Category', 'parameter_symbol': 'liq_category', 'parameter_value': liq_category},
                                    {'parameter_name': 'Liquid Name', 'parameter_symbol': 'liq_name', 'parameter_value': liq_name},
                                    {'parameter_name': 'Liquid Bulk Temp', 'parameter_symbol': 'tb', 'parameter_value': tb},
                                    {'parameter_name': 'Average Daily Liquid Surface Temp, R', 'parameter_symbol': 'tla', 'parameter_value': tla},
                                    {'parameter_name': 'Average Daily Min Liquid Surface Temp, R', 'parameter_symbol': 'tln', 'parameter_value': tln},
                                    {'parameter_name': 'Average Daily Max Liquid Surface Temp, R', 'parameter_symbol': 'tlx', 'parameter_value': tlx},
                                    {'parameter_name': 'Vapor Pressure at the Ave. Daily Min Liquid Surface Temp, R', 'parameter_symbol': 'pva', 'parameter_value': pva},
                                    {'parameter_name': 'Vapor Pressure at the Ave. Daily Max Liquid Surface Temp, R', 'parameter_symbol': 'pvn', 'parameter_value': pvn},
                                    {'parameter_name': 'Vapor Molecular Weight lb/lb-mole', 'parameter_symbol': 'pvx', 'parameter_value': pvx},
                                    {'parameter_name': 'RVP, psia', 'parameter_symbol': 'mv', 'parameter_value': mv},
                                    {'parameter_name': 'Grant', 'parameter_symbol': 'rvp', 'parameter_value': rvp},
                                    {'parameter_name': 'Constant A, unitless', 'parameter_symbol': 'const_a', 'parameter_value': const_a},
                                    {'parameter_name': 'Constant B, R', 'parameter_symbol': 'const_b', 'parameter_value': const_b}]
        }],

        'standing_losses': [{

            'value': ls,
            'parameter_estimates': [{'parameter_name': 'Standing Losses, lbs', 'parameter_symbol': 'ls', 'parameter_value': ls},
                                    {'parameter_name': 'Vapor Space Volume, ft^3', 'parameter_symbol': 'vv', 'parameter_value': vv},
                                    {'parameter_name': 'Vapor Density, lb/ft^3', 'parameter_symbol': 'wv', 'parameter_value': wv},
                                    {'parameter_name': 'Vapor Space Expansion Factor, 1/day', 'parameter_symbol': 'ke', 'parameter_value': ke},
                                    {'parameter_name': 'Vented Vapor Saturation Factor, unitless', 'parameter_symbol': 'ks', 'parameter_value': ks}]
        }],

        'working_losses': [{

            'value': lw,
            'parameter_estimates': [{'parameter_name': 'Working Losses, lbs/year', 'parameter_symbol': 'lw', 'parameter_value': lw},
                                    {'parameter_name': 'Net Working Loss Throughput, ft^3/year', 'parameter_symbol': 'vq', 'parameter_value': vq},
                                    {'parameter_name': 'Turnover Factor, unitless', 'parameter_symbol': 'kn', 'parameter_value': kn},
                                    {'parameter_name': 'Working Loss Product Factor, unitless', 'parameter_symbol': 'kp', 'parameter_value': kp},
                                    {'parameter_name': 'Vapor Density, lb/ft^3', 'parameter_symbol': 'wv', 'parameter_value': wv},
                                    {'parameter_name': 'Vent Setting Correction Factor, unitless', 'parameter_symbol': 'kb', 'parameter_value': kb}]
        }],

        'total_losses': [{

            'value': lt,
            'parameter_estimates': [{'parameter_name': 'working_loss, lbs', 'parameter_symbol': 'lw', 'parameter_value': lw},
                                    {'parameter_name': 'standing_loss, lbs', 'parameter_symbol': 'ls', 'parameter_value': ls},
                                    {'parameter_name': 'total_loss, lbs', 'parameter_symbol': 'lt', 'parameter_value': lt}]
        }],

        'vapor_space_volume': [{

            'value': vv,
            'parameter_estimates': [{'parameter_name': 'Vapor Space Volume', 'parameter_symbol': 'vv', 'parameter_value': vv},
                                    {'parameter_name': 'Tank Diameter', 'parameter_symbol': 'diameter', 'parameter_value': diameter},
                                    {'parameter_name': 'Vapor Space Outage', 'parameter_symbol': 'hvo', 'parameter_value': hvo}]
        }],

        'vapor_space_outage': [{

            'value': hvo,
            'parameter_estimates': [{'parameter_name': 'Vapor Space Outage', 'parameter_symbol': 'hvo', 'parameter_value': hvo},
                                    {'parameter_name': 'Tank Shell Height', 'parameter_symbol': 'hs', 'parameter_value': hs},
                                    {'parameter_name': 'Liquid Height', 'parameter_symbol': 'hl', 'parameter_value': hl},
                                    {'parameter_name': 'Roof Outage', 'parameter_symbol': 'hro', 'parameter_value': hro}]
        }],

        'roof_outage': [{

            'value': hro,
            'parameter_estimates': [{'parameter_name': 'Roof Outage', 'parameter_symbol': 'hro', 'parameter_value': hro},
                                    {'parameter_name': 'Tank Roof Height', 'parameter_symbol': 'hr', 'parameter_value': hr},
                                    {'parameter_name': 'Tank Shell Radius', 'parameter_symbol': 'rs', 'parameter_value': rs},
                                    {'parameter_name': 'Tank Cone Roof Slope', 'parameter_symbol': 'sr', 'parameter_value': sr}]
        }],

        'vapor_density': [{

            'value': wv,
            'parameter_estimates': [{'parameter_name': 'Vapor Density, lb/ft^3', 'parameter_symbol': 'wv', 'parameter_value': wv},
                                    {'parameter_name': 'Vapor Molecular Weight, lb/lb-mole', 'parameter_symbol': 'mv', 'parameter_value': mv},
                                    {'parameter_name': 'Vapor Pressure at Average Daily Liquid Surface Temp, psia', 'parameter_symbol': 'pva', 'parameter_value': pva},
                                    {'parameter_name': 'Ideal Gas Constant', 'parameter_symbol': 'r', 'parameter_value': r},
                                    {'parameter_name': 'Average Vapor Temperature', 'parameter_symbol': 'tv', 'parameter_value': tv},
                                    {'parameter_name': 'Tank Roof Surface Solar Absorptance', 'parameter_symbol': 'alpha_s', 'parameter_value': alpha_s},
                                    {'parameter_name': 'Tank Shell Surface Solar Absorptance', 'parameter_symbol': 'alpha_r', 'parameter_value': alpha_r},
                                    {'parameter_name': 'Average Daily Total Insolation Factor', 'parameter_symbol': 'i_value', 'parameter_value': i_value}]
        }],

        'vapor_space_expansion_factor': [{

            'value': ke,
            'parameter_estimates': [{'parameter_name': 'Vapor Space Expansion Factor, 1/day', 'parameter_symbol': 'ke', 'parameter_value': ke},
                                    {'parameter_name': 'Average Daily Vapor Presure Temp. Range, R', 'parameter_symbol': 'del_tv', 'parameter_value': del_tv},
                                    {'parameter_name': 'Average Daily Vapor Pressure Range, psi', 'parameter_symbol': 'delta_pv', 'parameter_value': delta_pv},
                                    {'parameter_name': 'Breather Vent Pressure Setting Range, psig', 'parameter_symbol': 'delta_pb', 'parameter_value': delta_pb},
                                    {'parameter_name': 'Atmospheric Pressure, psi', 'parameter_symbol': 'pa', 'parameter_value': pa},
                                    {'parameter_name': 'Vapor Pressure at Average Daily Liquid Surface Temp, psia', 'parameter_symbol': 'pva', 'parameter_value': pva},
                                    {'parameter_name': 'Average Daily Liquid Surface Temp, R', 'parameter_symbol': 'tla', 'parameter_value': tla}]
        }],

        'vented_vapor_space_saturation_factor': [{

            'value': ks,
            'parameter_estimates': [{'parameter_name': 'Vented Vapor Saturation Factor, unitless', 'parameter_symbol': 'ks', 'parameter_value': ks},
                                    {'parameter_name': 'Vapor Pressure at Average Daily Liquid Surface Temp, psia', 'parameter_symbol': 'pva', 'parameter_value': pva},
                                    {'parameter_name': 'Vapor Space Outage, ft', 'parameter_symbol': 'hvo', 'parameter_value': hvo}]
        }],
        
        'data_inputs': [{
            
            'value': None,
            'parameter_estimates': [{'parameter_name': 'hs_shellHeight', 'parameter_symbol': 'g', 'value': data[9]['data']['shellHeight']}, 
                                    {'parameter_name': 'tankLength', 'parameter_symbol': 'g', 'value': data[9]['data']['tankLength']},
                                    {'parameter_name': 'tankDiameter', 'parameter_symbol': 'g', 'value': data[9]['data']['tankDiameter']},
                                    {'parameter_name': 'rr_tankDomeRoofRadius', 'parameter_symbol': 'g', 'value': data[9]['data']['tankDomeRoofRadius']},
                                    {'parameter_name': 'isTankDomeRoofRadiusKnown', 'parameter_symbol': 'g', 'value': data[9]['data']['isTankDomeRoofRadiusKnown']},
                                    {'parameter_name': 'isTankOnVaporBalance', 'parameter_symbol': 'g', 'value': data[9]['data']['isTankOnVaporBalance']},
                                    {'parameter_name': 'roofType', 'parameter_symbol': 'g', 'value': data[9]['data']['roofType']},
                                    {'parameter_name': 'tankOrientation', 'parameter_symbol': 'g', 'value': data[9]['data']['tankOrientation']},
                                    {'parameter_name': 'throughput', 'parameter_symbol': 'g', 'value': data[9]['data']['throughput']},
                                    {'parameter_name': 'productType', 'parameter_symbol': 'g', 'value': data[9]['data']['productType']},
                                    {'parameter_name': 'breatherVentSetting', 'parameter_symbol': 'g', 'value': data[9]['data']['breatherVentSetting']},
                                    {'parameter_name': 'breatherPressureSetting', 'parameter_symbol': 'g', 'value': data[9]['data']['breatherPressureSetting']},
                                    {'parameter_name': 'breatherVaccumSetting', 'parameter_symbol': 'g', 'value': data[9]['data']['breatherVaccumSetting']},
                                    {'parameter_name': 'location', 'parameter_symbol': 'g', 'value': data[9]['data']['location']},
                                    {'parameter_name': 'shellShade', 'parameter_symbol': 'g', 'value': data[9]['data']['shellShade']},
                                    {'parameter_name': 'shellCondition', 'parameter_symbol': 'g', 'value': data[9]['data']['shellCondition']},
                                    {'parameter_name': 'roofShade', 'parameter_symbol': 'g', 'value': data[9]['data']['roofShade']},
                                    {'parameter_name': 'roofCondition', 'parameter_symbol': 'g', 'value': data[9]['data']['roofCondition']},
                                    {'parameter_name': 'isShellInsulated', 'parameter_symbol': 'g', 'value': data[9]['data']['isShellInsulated']},
                                    {'parameter_name': 'isRoofInsulated', 'parameter_symbol': 'g', 'value': data[9]['data']['isRoofInsulated']},
                                    {'parameter_name': 'productClassString', 'parameter_symbol': 'g', 'value': data[9]['data']['productClassString']},
                                    {'parameter_name': 'productNameString', 'parameter_symbol': 'g', 'value': data[9]['data']['productNameString']},
                                    {'parameter_name': 'roofSlope', 'parameter_symbol': 'g', 'value': data[9]['data']['roofSlope']},
                                    {'parameter_name': 'isAverageLiquidHeightKnown', 'parameter_symbol': 'g', 'value': data[9]['data']['isAverageLiquidHeightKnown']},
                                    {'parameter_name': 'averageLiquidHeight', 'parameter_symbol': 'g', 'value': data[9]['data']['averageLiquidHeight']}, 
                                    {'parameter_name': 'rvp_crudeOils', 'parameter_symbol': 'g', 'value': data[9]['data']['rvp_crudeOils']},
                                    {'parameter_name': 'tankName', 'parameter_symbol': 'g', 'value': data[9]['data']['tankName']}, 
                                    {'parameter_name': 'assetNumber', 'parameter_symbol': 'g', 'value': data[9]['data']['assetNumber']}]
        }]

    }
    
    return return_data
