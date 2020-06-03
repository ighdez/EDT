from edt_functions import *

# Example 2: A more complex D-efficient design

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
            'par':		[.01,.02]}

Att2 = {	'name':		'B',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[.02,.05]}

Att3 = {	'name':		'C',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[.05,.1]}

Att4 = {	'name':		'D',
            'levels':	[1,2,3,4],
            'coding':	1,
            'par':		[.0025]}

Att5 = {	'name':		'E',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[-.1,-.2]}

Att6 = {	'name':		'F',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[-.05,-.1]}

Att7 = {	'name':		'G',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[0,0]}

Att8 = {	'name':		'H',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[-.15,-.2]}

Att9 = {	'name':		'I',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[-.1,-.15]}

Att10 = {	'name':		'J',
            'levels':	[0,1,2],
            'coding':	2,
            'par':		[-.02,-.06]}

# Set conditions
cond = 	[	"if D>2, H>0 and I>0 and J>0",
            "if D<3, H<2 and I<2 and J<2"]

# Merge all attributes in a single list
Attlist = [Att1,Att2,Att3,Att4,Att5,Att6,Att7,Att8,Att9,Att10]

# Set number of alternatives per choice set and number of choice sets
nalternatives = 2
nchoicesit = 72

# Options
options = 	{	'optout': 1,
                'asc': [3],
                'ascpar': [-1],
                'seed': 666,
                # 'iterlim': 0,
                # 'noimprovlim': 0,
                'timelim': 0.1,
                'nblocks': 9,
                'cond': cond}

design = effdesign(Attlist,nalternatives,nchoicesit,options)

to_export = design['optimal.des']

to_export.to_excel("efficient_design_cond.xlsx",index=False)