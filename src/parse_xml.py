import xml.etree.ElementTree as ET 
from xml.dom import minidom
from bs4 import BeautifulSoup as bs
import yaml



def parse_xml():

    with open('../config.yaml', 'r') as file:
        prime_service = yaml.safe_load(file)

    print(prime_service)

    root = bs(prime_service['rla_xml_file'], 'xml')
    rla = root.new_tag('rla-config')
    root.append(rla)


    globalI = root.new_tag('globalindex')
    rla.append(globalI)

    g_val = root.new_tag('value')

    st = " "

    for i in range(len(prime_service['header'])):
        st += (str(i) + ' ' )

    g_val.string = st
    #g_val.children = [ i for i in range(len(prime_service['header']))]

    print(len(prime_service['input_files']))
    globalI.append(g_val)

    weights = root.new_tag('weights')
    rla.append(weights)


    w_val = root.new_tag('value')

    st_n = " "

    for i in range(len(prime_service['header'])):
        st_n += ('1' + ' ' )

    w_val.string = st_n

    weights.append(w_val)

    input_dict = prime_service['sample_output']

    i = 0

    for key in input_dict:
        i += 1
        filename = input_dict[key]
        dt = " "
        dt += ("dataset" + str(i))
        dataset = root.new_tag('dataset',attrs={"name":dt,"id":i})
        rla.append(dataset)

        d_val = root.new_tag('value')
        d_val.string = str(filename)
        dataset.append(d_val)

        dataset_index = root.new_tag("dataset_index")
        dataset_index_val = root.new_tag('value')
        dataset_index_val.string = st

        dataset_index.append(dataset_index_val)
        dataset.append(dataset_index)


    version_conf = root.new_tag("version-config-param",attrs={"id":"ComparisonGroup"})


    comp = prime_service['COMPARISON_ATTRIBUTE']
        
    for i in range(len(comp)):
        comparison = root.new_tag("comparison",attrs={"id":i})
        dist_method = root.new_tag("dist_calc_method")
        dist_val = root.new_tag("value")
        dist_val.string = "1"
        dist_method.append(dist_val)
        comparison.append(dist_method)
        attribute_indices = root.new_tag("comparing_attribute_indices",attrs={"id":"1"})
        attribute_indices_val = root.new_tag("value")
        attribute_indices_val.string = str(comp[i])
        attribute_indices.append(attribute_indices_val)
        comparison.append(attribute_indices)

        version_conf.append(comparison)

    thresh = root.new_tag("threshold")
    thresh_val = root.new_tag("value")
    thresh_val.string = str(prime_service['THRESHOLD'])
    thresh.append(thresh_val)

    version_conf.append(thresh)


    priority = root.new_tag("priority")
    pri_val = root.new_tag("value")
    pri_val.string = str(1)
    priority.append(pri_val)

    version_conf.append(priority)

    block = root.new_tag("block")
    index = root.new_tag("index")
    value_index = root.new_tag("value")
    value_index.string = str(prime_service['BLOCKING_ATTRIBUTE'])
    index.append(value_index)
    block.append(index)

    length = root.new_tag("length")
    value_len = root.new_tag("value")
    value_len.string = "5"
    length.append(value_len)
    block.append(length)

    type = root.new_tag("type")
    value_type = root.new_tag("value")
    value_type.string = "0"
    type.append(value_type)
    block.append(type)

    version_conf.append(block)
    rla.append(version_conf)

    output_function = root.new_tag("output_function",attrs={"id":"1"})
    output_filename = root.new_tag("output_filename")
    value_file = root.new_tag("value")
    value_file.string = "./RLA_CL_EXTRACT/output/output.txt"
    output_filename.append(value_file)
    output_function.append(output_filename)
    version_conf.append(output_function)

    results_log = root.new_tag("results_logging")
    filename_log = root.new_tag("filename")
    value_log = root.new_tag("value")
    value_log.string = "results.csv"
    filename_log.append(value_log)
    results_log.append(filename_log)

    version_conf.append(results_log)
        
    soup_string = str(root)
    tree = ET.XML(soup_string)

    from xml.dom import minidom

    xmlstr = minidom.parseString(ET.tostring(tree)).toprettyxml(indent="   ")
    with open("../data/config_test.xml", "w") as f:
        f.write(xmlstr)

    # with open ("../data/config_test.xml", "wb") as files :
    #     files.write()