import requests as r
import json

POST_JSON = {
    'shellHeight':20, 
    'tankDiameter':15, 
    'tankDomeRoofRadius':9, 
    'breatherPressureSetting':0.03,
    'breatherVaccumSetting':-0.03,
    'tankLength':10,
    'throughput':100000, 
    'roofSlope':0.0625, 
    'averageLiquidHeight':10, 
    'rvp_crudeOils':5,
    'productMW':50,
    'isTankDomeRoofRadiusKnown':False, 
    'isTankOnVaporBalance':False, 
    'isRoofInsulated':False,
    'isShellInsulated':False,   
    'isAverageLiquidHeightKnown':False,
    'isCustomProduct':False,
    'roofType':'Cone',
    'tankOrientation':'Vertical', 
    'productType':'Crude Oils', 
    'breatherVentSetting':'Default',
    'location':'Long Beach, CA',
    'shellShade':'White', 
    'shellCondition':'Average', 
    'roofShade':'White',
    'roofCondition':'Average', 
    'productClassString':'Crude Oils',
    'productNameString':'Crude oil (RVP 5)'}

url = 'https://tanks-409d-api.herokuapp.com/apiv001/vfrt'
response = r.post(url, json=POST_JSON)
data = response.text
data_json = json.loads(data)
data_json
