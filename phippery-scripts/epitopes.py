
#EPITOPES = {
#    "NTD" : {"limits": [288, 307], "binding threshold": 40},
#    "CTD" : {"limits": [545, 780], "binding threshold": 200},
#    "CTD-1" : {"limits": [545, 585], "binding threshold": 120},
#    "CTD-2" : {"limits": [620, 680], "binding threshold": 120},
#    "FP" : {"limits": [805, 835], "binding threshold": 100},
#    "HR2" : {"limits": [1135,1170], "binding threshold": 150}
#}

PEPTIDE_FLANK = 15
EPITOPES = {
    "NTD" : {
        "limits": [285, 305], 
        "binding threshold": 40,
        #"left border": 289-288,
        #"right border": 306-288
        "left border": 289,
        "right border": 306
    },
    "CTD" : {
        "limits": [545, 690], 
        "binding threshold": 200,
        "left border": 550,
        "right border": 685
    },
    "CTD-1" : {
        "limits": [545, 580], 
        "binding threshold": 30,
        "left border": 550,
        "right border": 573
    },
    #"CTD-2" : {
    #    "limits": [620, 650], 
    #    "binding threshold": 30,
    #    "left border": 625,
    #    "right border": 646
    #},
    #"CTD-3" : {
    #    "limits": [650, 680], 
    #    "binding threshold": 30,
    #    "left border": 659,
    #    "right border": 677
    #},
    "FP" : {
        "limits": [805, 835], 
        "binding threshold": 100,
        "left border": 810,
        "right border": 829
    },
    "HR2" : {
        "limits": [1135,1170], 
        "binding threshold": 150,
        "left border": 1141,
        "right border": 1163
    }
}
