import json
import math
from sqlite3 import connect
import sqlite3
import pandas as pd
import requests as r
from flask import Flask, jsonify, request
from verticalFixedRoofTank import VerticalFixedRoofTank
from vfrtRun import returnData1

app = Flask(__name__)

@app.route('/is_alive')
def is_alive():
    return 'This is Tanks 4.09_d on the Web'

@app.route('/apiv001/vfrt', methods=['POST'])
def main():

    POST_DATA = request.get_json()

    INPUT_shellHeight = POST_DATA['shellHeight']
    INPUT_tankDiameter = POST_DATA['tankDiameter'] 
    INPUT_tankDomeRoofRadius = POST_DATA['tankDomeRoofRadius']
    INPUT_breatherPressureSetting = POST_DATA['breatherPressureSetting']
    INPUT_breatherVaccumSetting = POST_DATA['breatherVaccumSetting']
    INPUT_tankLength = POST_DATA['tankLength']
    INPUT_throughput = POST_DATA['throughput']
    INPUT_roofSlope = POST_DATA['roofSlope']
    INPUT_averageLiquidHeight = POST_DATA['averageLiquidHeight']
    INPUT_rvp_crudeOils = POST_DATA['rvp_crudeOils']
    INPUT_productMW = POST_DATA['productMW']
    INPUT_isTankDomeRoofRadiusKnown = POST_DATA['isTankDomeRoofRadiusKnown'] 
    INPUT_isTankOnVaporBalance = POST_DATA['isTankOnVaporBalance']
    INPUT_isRoofInsulated = POST_DATA['isRoofInsulated']
    INPUT_isShellInsulated = POST_DATA['isShellInsulated']  
    INPUT_isAverageLiquidHeightKnown = POST_DATA['isAverageLiquidHeightKnown']
    INPUT_isCustomProduct = POST_DATA['isCustomProduct']
    INPUT_roofType = POST_DATA['roofType']
    INPUT_tankOrientation = POST_DATA['tankOrientation']
    INPUT_productType = POST_DATA['productType']
    INPUT_breatherVentSetting = POST_DATA['breatherVentSetting']
    INPUT_location = POST_DATA['location']
    INPUT_shellShade = POST_DATA['shellShade']
    INPUT_shellCondition = POST_DATA['shellCondition']
    INPUT_roofShade = POST_DATA['roofShade']
    INPUT_roofCondition = POST_DATA['roofCondition'] 
    INPUT_productClassString = POST_DATA['productClassString']
    INPUT_productNameString = POST_DATA['productNameString']
    INPUT_assetNumber = POST_DATA['assetNumber']
    INPUT_tankName = POST_DATA['tankName']

    obj = returnData1(

            shellHeight=INPUT_shellHeight, 
            tankDiameter=INPUT_tankDiameter, 
            tankDomeRoofRadius=INPUT_tankDomeRoofRadius, 
            breatherPressureSetting=INPUT_breatherPressureSetting,
            breatherVaccumSetting=INPUT_breatherVaccumSetting,
            tankLength=INPUT_tankLength,
            throughput=INPUT_throughput, 
            roofSlope=INPUT_roofSlope, 
            averageLiquidHeight=INPUT_averageLiquidHeight, 
            rvp_crudeOils=INPUT_rvp_crudeOils,
            productMW=INPUT_productMW,
            isTankDomeRoofRadiusKnown=INPUT_isTankDomeRoofRadiusKnown, 
            isTankOnVaporBalance=INPUT_isTankOnVaporBalance, 
            isRoofInsulated=INPUT_isRoofInsulated,
            isShellInsulated=INPUT_isShellInsulated,   
            isAverageLiquidHeightKnown=INPUT_isAverageLiquidHeightKnown,
            isCustomProduct=INPUT_isCustomProduct,
            roofType=INPUT_roofType,
            tankOrientation=INPUT_tankOrientation, 
            productType=INPUT_productType, 
            breatherVentSetting=INPUT_breatherVentSetting,
            location=INPUT_location,
            shellShade=INPUT_shellShade, 
            shellCondition=INPUT_shellCondition, 
            roofShade=INPUT_roofShade,
            roofCondition=INPUT_roofCondition, 
            productClassString=INPUT_productClassString,
            productNameString=INPUT_productNameString, 
            tankName=INPUT_tankName, 
            assetNumber=INPUT_assetNumber)

    return jsonify(obj)

if __name__ == "__main__":
    app.run()
