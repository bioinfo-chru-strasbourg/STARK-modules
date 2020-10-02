import yaml

#Load the yaml files
with open('STARK.docker-compose.yml') as fp:
    data = yaml.load(fp, Loader=yaml.FullLoader)
with open('STARK.DECoN.docker-compose.yml') as fp:
    data1 = yaml.load(fp, Loader=yaml.FullLoader)
#Add the resources from test2.yaml to test1.yaml resources
# for i in data1['resources']:
#     print i,data1['resources'][i]
#     data['resources'].update({i:data1['resources'][i]})
#create a new file with merged yaml
with open('/tmp/store_file.yaml', 'w') as file:
    yaml.dump(data,file)