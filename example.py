from edt_functions import *

# Example 1: A simple D-efficient designs

# Define attributes
# Each attribute is a dictionary that contains:
#
# Key			Type			Description
# 'name'		String			Attribute name
# 'levels'		List			Attribute levels
# 'coding'		Integer			Coding. 1 = Quantitative / 2 = Dummy
# 'par'			List			Prior parameters

Att1 = {	'name':		'A',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[-0.1,-0.2]}

Att2 = {	'name':		'B',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[0.1,0.15]}

Att3 = {	'name':		'C',
            'levels':	[0,1],
            'coding':	2,
            'par':		[-0.1]}

Att4 = {	'name':		'D',
            'levels':	[0,1],
            'coding':	2,
            'par':		[-0.4]}

Att5 = {	'name':		'E',
            'levels':	[0,1],
            'coding':	2,
            'par':		[-0.1]}

Att6 = {	'name':		'F',
            'levels':	[0,1],
            'coding':	2,
            'par':		[-0.4]}

Att7 = {	'name':		'G',
            'levels':	[100,250,500,750],
            'coding':	1,
            'par':		[-0.004]}

# Merge all attributes in a single list
Attlist = [Att1,Att2,Att3,Att4,Att5,Att6,Att7]

# Set number of alternatives per choice set and number of choice sets
nalternatives = 2
nchoicesit = 60

# Options
options = 	{	'optout': 1,
                'asc': [],
                'ascpar': [-1],
                'seed': 666,
                # 'iterlim': 0,
                # 'noimprovlim': 0,
                'timelim': 1,
                'nblocks': 10}

design = effdesign(Attlist,nalternatives,nchoicesit,options)

to_export = design['optimal.des']

to_export.to_excel("efficient_design.xlsx",index=False)