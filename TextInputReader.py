from configparser import ConfigParser
parser = ConfigParser()
parser.read('/Users/bradyhislop/Box/ClusteringGUI_python/config.ini')

sections = parser.sections()

for i in range(len(sections)):
	print(sections[i])


