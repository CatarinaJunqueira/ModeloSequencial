# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 10:12:24 2018

@author: catar
"""

import requests

def email_alert(Y, z, out_string):
    report = {}
    report["value1"] = Y
    report["value2"] = z
    out_string = out_string.replace("\n","<br>" )
    report["value3"] = out_string

    requests.post("https://maker.ifttt.com/trigger/experimento_modelointegrado_concluido/with/key/ho6u3RW6TW1B1DqQe1g-V4zwuX86E_Hgjm8lvUOmy2f/", data=report)    

    return