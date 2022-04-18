
import csv, sys

ran = str(sys.argv[1])
serv = "phagetb_"
dir = serv+ran

csvFile = open('/usr1/webserver/cgidocs/tmp/phagetb/phagetb_'+ran+'/outputfile_' +ran+'.csv')#enter the csv filename
with open('/usr1/webserver/cgidocs/tmp/phagetb/phagetb_'+ran+'/output.html'+ran, 'w') as html: #enter the output filename
    with open('/usr1/webserver/cgidocs/tmp/phagetb/phagetb_'+ran+'/outputfile_' +ran+'.csv') as file:
     r = 0
     for row in file.read().split('\n')[:-1]:
        if r == 0:
            html.write('\t<thead class=\"thead-dark\">\r\t\t<tr>\r')
            html.write('\t\t\t<th data-sortable="true">' + "Query No." + '</th>\r')
            for col in row.split(','):
                html.write('\t\t\t<th data-sortable="true">' + col + '</th>\r')
            html.write('\t\t</tr>\r\t</thead>\r')
        elif r == 1:
            html.write('\t\t<tbody>\r')
            html.write('\t\t<tr>\r')
            html.write('\t\t\t<td>' + str(r) + '</td>\r')
            for col in row.split(','):
                html.write('\t\t\t<td>' + col + '</td>\r')
            html.write('\t\t</tr>\r')
        else:
            html.write('\t\t<tr>\r')
            html.write('\t\t\t<td>' + str(r) + '</td>\r')
            for col in row.split(','):
                html.write('\t\t\t<td>' + col + '</td>\r')
            html.write('\t\t</tr>\r')
        r += 1
    html.write('\t\t</tbody>\r')

