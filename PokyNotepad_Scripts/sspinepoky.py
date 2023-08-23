# -----------------------------------------------------------------------------
# ssPINE-POKY dialog
# -----------------------------------------------------------------------------
#
# Developed by Woonghee Lee, Andrea Estefania Lopez Giraldo, Mehdi Rahimi
# e-mail: woonghee.lee@ucdenver.edu
# Department of Chemistry, University of Colorado Denver
#
# Last updated: July 6, 2023
#
# ssPINE user group: https://groups.google.com/g/pinenmr-user-group
#
# -----------------------------------------------------------------------------
# CITATION:
#
#    ssPINE-POKY TBA
#
#    ssPINE: Probabilistic Algorithm for Automated Chemical Shift Assignment 
#    of Solid-State NMR Data from Complex Protein Systems. 
#    Adilakshmi Dwarasala, Mehdi Rahimi, John L Markley, Woonghee Lee. (2022) 
#    Membranes. 12(9): 834.
#    https://doi.org/10.3390/membranes12090834
# -----------------------------------------------------------------------------
#
# WHAT TO ADD? check if resonances without peaks and ask.
#

import os
import tkinter

import poky
import sputil
import tkutil
import myseq
import pokynmr
import zipfile
from sys import platform

# for token
import time

import tempfile
import requests
import serverside

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

sspine_types = ("2D-NCA", "2D-NCO", "2D-CC", "2D-NCACB", "NCACB",
              "NCACX", "NCOCACB", "NCOCX", "CAN(CO)CX", "NCACO",
              "CANCO", "NCOCA", "CAN(CO)CA", "CAN(CO)CACB")

expected_counts = (1, 1, 6, 1, 2, 3, 4, 1, 3, 1, 2, 2, 3, 2)
sspine_field = ('file01', 'file15', 'file02', 'file08', 'file03', 'file04',
              'file05', 'file06', 'file07', 'file09', 'file10', 'file11',
              'file12', 'file13')
sspine_field2 = ('exp1', 'exp15', 'exp2', 'exp8', 'exp3', 'exp4', 'exp5', 'exp6',
               'exp7', 'exp9', 'exp10', 'exp11', 'exp12', 'exp13')

def getFirstKey(item):
    return item[0]


class sspine_dialog(tkutil.Dialog, tkutil.Stoppable):

    def __init__(self, session):

        self.session = session
        self.verbose = 1
        self.msg = ''

        tkutil.Dialog.__init__(self, session.tk, 'ssPINE-POKY automation')

        # User name
        user_name_ef = tkutil.entry_field(self.top, 'Your name: ', initial='',
                                          width=40)
        user_name_ef.frame.pack(side='top', anchor='w', expand=0)
        self.user_name = user_name_ef.variable

        # User email
        user_email_ef = tkutil.entry_field(self.top, 'Your email: ', initial='',
                                           width=40)
        user_email_ef.frame.pack(side='top', anchor='w')
        self.user_email = user_email_ef.variable
        sep = tkinter.Frame(self.top, height=2, bd=1, relief="ridge")
        sep.pack(fill="both", padx=5, pady=5, side='top', expand=0)

        # Sequence file
        user_seq_ff = tkutil.file_field2(self.top, 'Sequence file: ', 'sspinenmr',
                                         '3-letter-code sequence file with indices (*.seq);; ' +
                                         '1-letter-code sequence file (*.fasta);; ' +
                                         '3-letter-code sequence file without indices (*.txt)',
                                         default_ext='.seq', width=35, session=session)
        user_seq_ff.frame.pack(side='top', anchor='w')
        self.user_seq = user_seq_ff.variable
        sep = tkinter.Frame(self.top, height=2, bd=1, relief="ridge")
        sep.pack(fill="both", padx=5, pady=5, side='top', expand=0)
        # check if file exists
        if session.project != None:
            if session.project.poky_directory != '' and \
                    os.path.exists(session.project.save_path):
                seq_path = session.project.save_path + '.seq'
                if os.path.exists(seq_path):
                    user_seq_ff.variable.set(seq_path)

        # Scrolling List
        self.listbox = tkutil.scrolling_list(self.top,
                                             'Spectra to run ssPINE / ssPINE-POKY', 5)
        self.listbox.frame.pack(side='top', fill='both', expand=1)

        sep = tkinter.Frame(self.top, height=2, bd=1, relief="ridge")
        sep.pack(fill="both", padx=5, pady=5, side='top', expand=0)

        # spectrum and type
        spec_frame = tkinter.Frame(self.top, width=600, height=25)
        spec_frame.pack(side='top', anchor='w')

        self.sp = sputil.spectrum_combo_menu(
            session, spec_frame, 'Select spectrum: ')
        self.sp.frame.pack(side='top', anchor='nw', expand=0)
        self.selected_type = tkutil.combo_menu(spec_frame, "Experiment type",
                                               sspine_types, '2D-NCA', visible_rows=15)
        self.selected_type.frame.pack(side='top', anchor='nw')

        sep = tkinter.Frame(self.top, height=2, bd=1, relief="ridge")
        sep.pack(fill="both", padx=5, pady=5, side='top', expand=0)

        # buttons
        self.br = tkutil.button_row(self.top,
                                    ('Add', self.add_cb),
                                    ('Delete', self.delete_cb),
                                    ('Clear', self.clear_cb),
                                    ('Submit', self.submit_cb),
                                    ('ssPINE Web', self.web_cb),
                                    ('User group', self.user_group_cb),
                                    ('Close', self.close_cb),
                                    )
        self.br.frame.pack(side='top', fill='both', expand=0)
        tkutil.create_hint(self.br.buttons[0],
                           'Add a spectrum to run ssPINE or ssPINE-POKY importer. ' +
                           'Experiment type should match with your spectrum.')
        tkutil.create_hint(self.br.buttons[1],
                           'Delete selected spectrum in the list.')
        tkutil.create_hint(self.br.buttons[2],
                           'Clear all spectra in the list.')
        tkutil.create_hint(self.br.buttons[3],
                           'Run a ssPINE automated assignment without visiting the ' +
                           'ssPINE web page.')
        tkutil.create_hint(self.br.buttons[4],
                           'Export peak lists and open ssPINE web page for manual submission.')
        tkutil.create_hint(self.br.buttons[5],
                           'Visit ssPINE user group to resolve issues.')

        sep = tkinter.Frame(self.top, height=2, bd=1, relief="ridge")
        sep.pack(fill="both", padx=5, pady=5, side='top', expand=0)

        #
        self.ps_label_text = tkinter.StringVar()
        self.ps_label_text.set("- ssPINE-POKY: ssPINE-POKY importer -")
        self.ps_label = tkinter.Label(self.top, textvariable=self.ps_label_text,
                                      anchor='w', justify='left')
        self.ps_label.pack(side='top', anchor='w')

        # token
        # User name
        token_frame = tkinter.Frame(self.top, width=600, height=25)
        token_frame.pack(side='top', anchor='w')
        token_ef = tkutil.entry_field(token_frame, 'Key: ', initial='', width=27)
        token_ef.frame.pack(side='left', anchor='w', expand=0)
        self.token = token_ef.variable
        self.br2 = tkutil.button_row(token_frame,
                                     ('Check', self.check_wrap_cb),
                                     ('Browse...', self.sspine_browse_cb),
                                     ('Web...', self.web_report_cb),
                                     )
        self.br2.frame.pack(side='top', anchor='e', expand=0)
        tkutil.create_hint(self.br2.buttons[0],
                           'Check ssPINE job status and import the results.')
        tkutil.create_hint(self.br2.buttons[1],
                           'Import the results in the local machine.')
        tkutil.create_hint(self.br2.buttons[2],
                            'Open ssPINE web report page.')

        # tolerances
        tol_frame = tkinter.Frame(self.top)
        tol_frame.pack(side='top', fill='both')

        self.ntol = tkutil.entry_field(tol_frame, 'N tol.:',
                                       initial=0.500, width=5)
        self.ntol.frame.grid(row=0, column=0)
        self.ntol_ppm = self.ntol.variable

        self.ctol = tkutil.entry_field(tol_frame, ', C tol.:',
                                       initial=0.500, width=5)
        self.ctol.frame.grid(row=0, column=2)
        self.ctol_ppm = self.ctol.variable

        # self.htol = tkutil.entry_field(tol_frame, ', H tol.:',
        #                                initial=0.050, width=5)
        # #self.htol.frame.grid(row=0, column=4)
        # self.htol_ppm = self.htol.variable
        
        #
        sep = tkinter.Frame(self.top, height=2, bd=1, relief="ridge")
        sep.pack(fill="both", padx=5, pady=5, side='top', expand=0)
        self.notice_text = tkinter.StringVar()
        self.notice_text.set("* Check peak counts if ssPINE " +
                             "keeps failing (Use 'pv' tool).")
        self.notice = tkinter.Label(self.top, textvariable=self.notice_text,
                                    anchor='w', justify='left')
        self.notice.pack(side='top', anchor='w')

    # -----------------------------------------------------------------------------
    #
    # sspine_params = [['',''],['',''],['',''],['',?],]
    def auto_setup(self, sspine_params):
        # 'user_name', 'user_email'
        # 'pre_assign', 'selective_label', 'cs_rosetta', 'seq_file',
        # 'HSQC (N15)', ....
        self.verbose = 0
        self.msg = ''
        self.listbox.clear()
        for key, value in sspine_params:
            if key == 'user_name':
                self.user_name.set(value)
            elif key == 'user_email':
                self.user_email.set(value)
            elif key == 'pre_assign':
                self.preassign.set(value)
            elif key == 'selective_label':
                self.labeling.set(value)
            elif key == 'cs_rosetta':
                self.csrosetta.set(value)
            elif key == 'seq_file':
                self.user_seq.set(value)
            else:  # experiments
                try:
                    #exp_no = sspine_types.index(key)
                    self.listbox.append('[%s]:%s' %
                                        (key, value.name), (key, value))
                except:
                    self.msg = 'Unknown experiment: ' + key
                    return
        self.submit_cb()
        if self.msg == '':
            self.msg = 'Submission succeeded.'
        self.verbose = 1
    # -----------------------------------------------------------------------------
    
    def add_cb(self):
        spec = self.sp.spectrum()
        spec_type = self.selected_type.get()
        # already in the list?
        iFound = 0

        for data in self.listbox.line_data:
            if data[0] == spec_type and data[1] == spec:
                iFound = 3
                break
            if data[0] == spec_type:
                iFound = 1
                break
            if data[1] == spec:
                iFound = 2
                break
        if iFound == 3 and spec_type == 'Selective_labeling':
            self.msg = '[%s]:%s already added.' % (spec_type, spec.name)
            if self.verbose == 1:
                self.session.show_message('Already added', self.msg)
        if iFound == 1 and spec_type != 'Selective_labeling':
            self.msg = '%s already added.' % (spec_type)
            if self.verbose == 1:
                self.session.show_message('Already added', self.msg)
            return
        if iFound == 2 and spec_type != 'Selective_labeling':
            self.msg = '%s already added.' % (spec.name)
            if self.verbose == 1:
                self.session.show_message('Already added', self.msg)
            return
        self.listbox.append('[%s]:%s' %
                            (spec_type, spec.name), (spec_type, spec))

    # -----------------------------------------------------------------------------
    #
    def delete_cb(self):
        self.listbox.delete_selected()
    # -----------------------------------------------------------------------------
    #

    def clear_cb(self):
        self.listbox.clear()
    # -----------------------------------------------------------------------------
    #

    def send_to_serverside(self, program, fields, contents):
        base = tempfile.gettempdir()
        files = []
        parameters = {}

        for item in fields:
            parameters[item[0]] = item[1]

        data_type_path = os.path.join(base, 'data_type.txt')
        if os.path.isfile(data_type_path):
            os.remove(data_type_path)
        files.append(data_type_path)

        for item in contents:
            file_number = item[0]
            try:
                type_index = sspine_field.index(file_number)
                exp_number = sspine_field2[type_index]
            except:
                if file_number == 'seqfile':
                    exp_number = 'seq'
                elif file_number == 'pdb':
                    exp_number = 'pdb_file'
                elif file_number == 'labeling':
                    exp_number = 'label_file'
                elif file_number == 'pre_assign':
                    exp_number = 'assign_file'

            file_type = parameters[exp_number]
            file_name = serverside.get_clean_filename(item[1])
            file_content = item[2]
            file_path = os.path.join(base, file_name)
            files.append(file_path)

            with open(file_path, 'w') as f:
                f.write(file_content)

            with open(data_type_path, 'a') as f:
                f.write(file_type + ',' + file_name + '\n')

        serverside_result = serverside.send(program, files, parameters)
        return serverside.get_jobid(serverside_result)

    # -----------------------------------------------------------------------------
    #
    def submit_cb(self):
        # sequence
        if not os.path.exists(self.user_seq.get()):
            self.msg = 'Please set your sequence file.'
            if self.verbose == 1:
                self.session.show_message('Sequence', self.msg)
            return

        readSeq = myseq.ReadSequenceFromFile(self.user_seq.get())
        oneline = ''
        sspine_seq = ''

        iSeqCount = len(readSeq)
        for seq in readSeq:
            if seq[1] == 'X':
                self.msg = 'Wrong sequence format: '
                if self.verbose == 1:
                    self.session.show_message('Sequence', self.msg + oneline)
                return
            oneline = oneline+seq[1]
            sspine_seq = sspine_seq + seq[1] + '\n'

        if iSeqCount < 5:
            self.msg = 'Your sequence length is only %d. ' + \
                'Please use different sequence file if it is incorrect.' + \
                str(iSeqCount)
            if self.verbose == 1:
                self.session.show_message('Sequence', self.msg)

        # check if available name and email exist
        szUserName = self.user_name.get()
        szUserEmail = self.user_email.get()

        if szUserName == '':
            szUserName = 'anonymous'

        if (szUserEmail == '') or (szUserEmail.find('@') == -1):
            szUserEmail = 'poky.nmr@gmail.com'

        fields = [['username', szUserName],
                  ['useremail', szUserEmail]]

        # make up spectra
        contents = [['seqfile', 'poky_seqfile.txt', sspine_seq]]
        fields.append(['seq', 'sequence'])

        # selective labeling should have selective labeling data
        iFound = 0
        szWrite = ''

        for (spec_type, spec) in self.listbox.line_data:
            try:
                iIdx = sspine_types.index(spec_type)
            except:
                self.msg = 'What is this experiment??? - ' + spec_type
                if self.verbose == 1:
                    self.session.show_message('Wrong experiment', self.msg)
                return

            peaklist = '      Assignment    '
            szHead = '?'
            for i in range(spec.dimension-1):
                peaklist = peaklist + '   w%d   ' % (i+1)
                szHead = szHead + '-?'
            peaklist = peaklist + '   w%d   ' % (spec.dimension)
            peaklist = peaklist + 'Data Height\n\n'

            if len(spec.peak_list()) == 0:
                self.msg = '%s does not have any peak.' % spec.name
                if self.verbose == 1:
                    self.session.show_message('No peaks', self.msg)
                return

            iCurrentCount = len(spec.peak_list())
            iExpectCount = expected_counts[iIdx]
            if iCurrentCount < int(float(iExpectCount) * 0.75):
                self.msg = '%s may need more peaks. Keep going anyway.' % spec.name
                if self.verbose == 1:
                    self.session.show_message('Warning', self.msg)

            peaklist2 = []
            for peak in spec.peak_list():
                if peak.data_height == None:
                    continue
                szWrite = '%10s' % szHead
                for i in range(spec.dimension):
                    szWrite = szWrite + ' %8.3f' % (peak.frequency[i])
                szWrite = szWrite + ' %16.3f' % (peak.data_height)
                szWrite = szWrite + '\n'
                peaklist2.append(szWrite)
            filtered_list = list(set(peaklist2))
            for peakitem in filtered_list:
                peaklist = peaklist + peakitem

            contents.append([sspine_field[iIdx], spec.name + '.list', peaklist])
            fields.append([sspine_field2[iIdx], sspine_types[iIdx]])

        token = self.send_to_serverside('ssPINE', fields, contents)
        self.token.set(token)

        if self.verbose == 1:
            self.session.show_message('ssPINE job sumitted.',
                                      'ssPINE job submitted. Check your result in this window ' +
                                      'or in your email inbox after a few minutes. \nPlease cite ' +
                                      'our ssPINE-POKY and ssPINE web server papers to support ' +
                                      'the maintenance (check with "n0").')
    # -----------------------------------------------------------------------------
    #
    def check_wrap_cb(self):
        self.check_cb(1)
    # -----------------------------------------------------------------------------
    #
    def sspine_browse_cb(self):
        self.check_cb(2)
    # -----------------------------------------------------------------------------
    #
    def web_report_cb(self):
        if self.token.get() == '':
            self.session.show_message('Error', 'ssPINE job ID not specified.')
            return
        url = 'https://poky.clas.ucdenver.edu/serverside/results/' + self.token.get() + \
            '/web_results/'
        import webbrowser
        webbrowser.open(url)

    # -----------------------------------------------------------------------------
    #
    def check_cb(self, clicked):
        if clicked == 2: # directory mode
            pine_dir = tkutil.load_directory_qt(self.session,
                                                'Select a ssPINE output directory',
                                                initdir=self.session.project.poky_directory)
            
            if pine_dir == '':
                return
            self.token.set(pine_dir)

        if self.token.get() == '':
            if clicked == 1:
                self.session.show_message('ssPINE job submitted.',
                                          'No submitted ssPINE job found.')
            return
        sspine_job_id = self.token.get()

        # retrieve ssPINE URL + token path.
        try:
            import urllib.request
            import urllib.error
            import urllib.parse
        except:
            self.session.show_message('urllib2',
                                      'You do not have openssl installed or set properly.\n' +
                                      'Please check your email instead.')
            return

        sspine_root_url = 'https://poky.clas.ucdenver.edu/serverside/results/%s/' % \
                        (sspine_job_id)
        sspine_result_url = sspine_root_url + 'results/'
        if clicked == 2:
             sspine_root_url = 'file://%s/' % (self.token.get())
             sspine_result_url = sspine_root_url

        sspine_package_url = sspine_result_url + 'all_files.zip'

        self.notice_text.set('Waiting ...')
        if clicked == 1:
            server_check = requests.head(sspine_package_url)
            if server_check.status_code != 200:
                self.session.show_message('Retrieve error',
                                          'It seems your job is not finished yet.')
                self.notice_text.set("* Check peak counts if ssPINE keeps" +
                                     " failing (Use 'pv' tool).")
                return
            else:
                self.session.show_message('ssPINE finished',
                                          'Your job is finished.')

        #if clicked == 2:
        sspine_ndpplot_url = sspine_result_url + 'ndpplot_ssPINE.ini'
        pecan_ndpplot_url = sspine_result_url + 'ndpplot_PECAN.ini'
        lacs_co_ndpplot_url = sspine_result_url + 'ndpplot_LACS_CACB_CO.ini'
        lacs_ca_ndpplot_url = sspine_result_url + 'ndpplot_LACS_CACB_CA.ini'
        lacs_cb_ndpplot_url = sspine_result_url + 'ndpplot_LACS_CACB_CB.ini'
        rci_s2_ndpplot_url = sspine_result_url + 'talosn_results/ndpplot_talosn_s2.ini'
        ini_path = os.path.join(os.path.expanduser("~"), 'ndpplot.ini')

        if self.session.show_message_yes_no('ssPINE finished',
                                            'Do you want to download and unpack ssPINE package?'):
            sspine_dir = os.path.join(self.session.project.poky_directory, 'ssPINE')
            sspine_job_dir = os.path.join(sspine_dir, sspine_job_id)
            sspine_path = os.path.join(sspine_job_dir, 'all_files.zip')

            try:
                if not os.path.exists(sspine_dir):
                    os.makedirs(sspine_dir)
                if not os.path.exists(sspine_job_dir):
                    os.makedirs(sspine_job_dir)

                sspine_response = urllib.request.urlopen(sspine_package_url, timeout=45)
                sspine_html = sspine_response.read()
                f = open(sspine_path, 'wb')
                f.write(sspine_html)
                f.close()

                zip_ref = zipfile.ZipFile(sspine_path, 'r')
                zip_ref.extractall(sspine_job_dir)
                zip_ref.close()

                self.notice_text.set('Results are downloaded')
            except:
                print(sspine_package_url)
                self.session.show_message("Error",
                                            "ssPINE results could not be retrieved or could " +
                                            "not be unpacked. Check ssPINE directory if you " +
                                            "have all_files.zip file.")
                self.notice_text.set('ssPINE results could not be retrieved')

        if self.session.show_message_yes_no('ssPINE finished',
                                            'Do you want to plot PECAN secondary structure prediction results?'):
            try:
                pecan_ndpplot = urllib.request.urlopen(pecan_ndpplot_url, timeout=10)
                temp = pecan_ndpplot.read().decode('utf-8')
                lines = temp.splitlines()
                # write ndpplot.ini
                f = open(ini_path, 'w')
                for line in lines:
                    f.write(temp+'\n')
                f.close()
                # open ndpplot
                pokynmr.run_npdplot(self.session)
            except:
                self.session.show_message('Retrieve error',
                                            'PECAN NDPPLOT results could not be retrieved.')
                pass

        if self.session.show_message_yes_no('ssPINE finished',
                                            'Do you want to plot ssPINE chemical shift assignment results?'):
            try:
                sspine_ndpplot = urllib.request.urlopen(sspine_ndpplot_url, timeout=10)
                temp = sspine_ndpplot.read().decode('utf-8')
                lines = temp.splitlines()
                if len(lines) > 5:
                    if str(lines[0]).find('[GENERAL_START]') != -1:
                        # write ndpplot.ini
                        f = open(ini_path, 'w')
                        for line in lines:
                            f.write(line+'\n')
                        f.close()
                        # open ndpplot
                        pokynmr.run_npdplot(self.session)
            except:
                self.session.show_message('Retrieve error',
                                        'ssPINE NDPPLOT results could not be retrieved.')
                pass

        if self.session.show_message_yes_no('ssPINE finished',
                                            'Do you want to plot LACS chemical shift validation results?'):
            try:
                lacs_co_ndpplot = urllib.request.urlopen(lacs_co_ndpplot_url, timeout=10)
                temp = lacs_co_ndpplot.read().decode('utf-8')
                lines = temp.splitlines()
                if len(lines) > 5:
                    if str(lines[0]).find('[GENERAL_START]') != -1:
                        # write ndpplot.ini
                        f = open(ini_path, 'w')
                        for line in lines:
                            f.write(line+'\n')
                        f.close()
                        # open ndpplot
                        pokynmr.run_npdplot(self.session)
                        time.sleep(3)
            except:
                pass

            try:
                lacs_cb_ndpplot = urllib.request.urlopen(lacs_cb_ndpplot_url, timeout=10)
                temp = lacs_cb_ndpplot.read().decode('utf-8')
                lines = temp.splitlines()
                if len(lines) > 5:
                    if str(lines[0]).find('[GENERAL_START]') != -1:
                        # write ndpplot.ini
                        f = open(ini_path, 'w')
                        for line in lines:
                            f.write(line+'\n')
                        f.close()
                        # open ndpplot
                        pokynmr.run_npdplot(self.session)
                        time.sleep(3)
            except:
                pass

            try:
                lacs_ca_ndpplot = urllib.request.urlopen(lacs_ca_ndpplot_url, timeout=10)
                temp = lacs_ca_ndpplot.read().decode('utf-8')
                lines = temp.splitlines()
                if len(lines) > 5:
                    if str(lines[0]).find('[GENERAL_START]') != -1:
                        # write ndpplot.ini
                        f = open(ini_path, 'w')
                        for line in lines:
                            f.write(line+'\n')
                        f.close()
                        # open ndpplot
                        pokynmr.run_npdplot(self.session)
            except:
                pass

        if self.session.show_message_yes_no('ss-PINE finished', 'Do you want to visualize ordered regions (RCI-S2 prediction)?'):
            try:
                rci_s2_ndpplot = urllib.request.urlopen(
                    rci_s2_ndpplot_url, timeout=10)
                temp = rci_s2_ndpplot.read().decode('utf-8')
                lines = temp.splitlines()
                if len(lines) > 5:
                    if str(lines[0]).find('[GENERAL_START]') == -1:
                        self.session.show_message(
                            'TALOS-N not finished', 'It seems TALOS-N job is not finished. Try again a few minutes later.')
                    else:
                        # write ndpplot.ini
                        f = open(ini_path, 'w')
                        for line in lines:
                            f.write(line+'\n')
                        f.close()
                        # open ndpplot
                        pokynmr.run_npdplot(self.session)
                        time.sleep(3)
                else:
                    self.session.show_message(
                        'TALOS-N not finished', 'It seems TALOS-N job is not finished. Try again a few minutes later.')
            except:
                self.session.show_message(
                    'TALOS-N not finished', 'It seems TALOS-N job is not finished. Try again a few minutes later.')
                pass

        
        if self.session.show_message_yes_no('ssPINE finished',
                                            'Do you want to create ssPINE labels by ssPINE-POKY ?'):
            # read - backbone, sidechain
            sspine_bb_url = sspine_result_url + 'protein_backbone_assignment.txt'
            sspine_sc_url = sspine_result_url + 'sidechain_table.txt'

            try:
                sspine_bb = urllib.request.urlopen(sspine_bb_url, timeout=10)
                sspine_bb_lines = sspine_bb.read().decode('utf-8').splitlines()  #It seems it is working [[1       MET    1.000  122.44  171.19   54.50   32.20   0.000    0.00    0.00   0.000    0.00    0.00   0.000    0.00    0.00   0.000 ]]
                sspine_sc = urllib.request.urlopen(sspine_sc_url, timeout=10)
                sspine_sc_lines = sspine_sc.read().decode('utf-8').splitlines()  #working also
            except:
                self.session.show_message(
                    'Retrieve error', 'ssPINE results could not be retrieved.')
                pass
            
            #print (sspine_bb_lines)

            if self.create_sspine_labels(sspine_bb_lines, sspine_sc_lines):  
                self.session.show_message(
                    'ssPINE-POKY', 'ssPINE labels created.')

                if self.session.show_message_yes_no('ssPINE finished',
                                                    'Do you want to assign peaks with ssPINE labels higher than 0.5 probability (A few seconds will be needed.) ?'):
                    self.assign_the_best() 

                self.session.show_message(
                    'ssPINE-POKY', 'Use ssPINE-POKY extensions to use ssPINE labels. Start with two-letter-code "pp". Floating ssPINE labels can be selected with two-letter-code "se". ')
            else:
                self.session.show_message(
                    'ssPINE label error', 'ssPINE labels could not be generated.')

        self.notice_text.set(
            "* Check peak counts if ssPINE keeps failing (Use 'pv' tool).")

    # -----------------------------------------------------------------------------

    def check_pattern(self, spec_type, w1, w2, w3, w4):  #check if the nucleus are from the same aa or not
        [shift, resonance, seq, seqidx, atom, prob] = w1
        [shift2, resonance2, seq2, seqidx2, atom2, prob2] = w2
        if w3 != None:
            [shift3, resonance3, seq3, seqidx3, atom3, prob3] = w3

        # sort by shift!
        if spec_type.find('4D') != -1:
            sort_list = sorted([w1, w2, w3, w4], key=getFirstKey)
            [w1, w2, w3, w4] = sort_list
        if spec_type.find('2D') != -1:
            sort_list = sorted([w1, w2], key=getFirstKey)
            [w1, w2] = sort_list
        else:
            sort_list = sorted([w1, w2, w3], key=getFirstKey)
            [w1, w2, w3] = sort_list
            
        [shift, resonance, seq, seqidx, atom, prob] = w1
        [shift2, resonance2, seq2, seqidx2, atom2, prob2] = w2
        if spec_type.find('2D') == -1:
            [shift3, resonance3, seq3, seqidx3, atom3, prob3] = w3
        if spec_type.find('4D') != -1:
            [shift4, resonance4, seq4, seqidx4, atom4, prob4] = w4

           # C-C  CX/O(i)-CX/O(i)  
        if (spec_type == '2D-CC'):
           if seqidx != seqidx2: #or abs(seqidx-seqidx2) > 1 :
               print(f'{spec_type}: {w1} {w2} {seqidx} {seqidx2} {atom} {atom2} 0')
               return 0
           #if atom not in ['CA', 'CB', 'CG', 'CG1','CG2', 'CD', 'CD1', 'CD2', 'CE', 'C'] or atom2 not in ['CA', 'CB', 'CG', 'CG1','CG2', 'CD', 'CD1', 'CD2', 'CE', 'C']:
           if atom[0] != 'C' or atom2[0] != 'C':
               print(f'{spec_type}: {w1} {w2} {seqidx} {seqidx2} {atom} {atom2} 0')
               return 0
           print(f'{spec_type}: {w1} {w2} {seqidx} {seqidx2} {atom} {atom2} 1')
           return 1

        # N - CX  Done
        if (spec_type == '2D-NCA'):  # N - CX It seems good
            if seqidx != seqidx2:
                return 0
            if (atom != 'CA' or atom2 != 'N') and (atom != 'N' or atom2 != 'CA'):
                return 0
            return 1
        
        #Done
        if spec_type == '2D-NCACB':  
            if seqidx != seqidx2:  #(if both are the same it is because the nuclei belong to the same aminoacid)
                return 0
            if (atom != 'N' or atom2 not in ['CA', 'CB']) and (atom2 != 'N' or atom not in ['CA', 'CB']):
                return 0
            return 1
        
        # N - C   Done
        if (spec_type == '2D-NCO'): # N - C
            if abs(seqidx - seqidx2) != 1:
                return 0
            if (atom != 'N' or atom2 != 'C') and (atom2 != 'N' or atom != 'C'):
                return 0
            return 1

        # N - CX - CX    Done
        #     if seqidx != seqidx3 or abs(seqidx-seqidx2) > 1 or seqidx < seqidx2:
        if spec_type == 'NCACB':
            if seqidx != seqidx2 or seqidx != seqidx3:
                return 0
            if sorted([atom, atom2, atom3]) != sorted(['N', 'CA', 'CB']):
                return 0
            return 1
        
        # N(i)-CA(i)-CX(i)  Done
        if spec_type == 'NCACX':
            if seqidx != seqidx2 or seqidx != seqidx3 or seqidx2 != seqidx3:
                #print(f'{spec_type}: {w1} {w2} {w3} {seqidx} {seqidx2} {seqidx3} {atom} {atom2} {atom3} 0')
                return 0
            try:
                n_idx = [atom, atom2, atom3].index('N')
            except:
                return 0
            try:
                ca_idx = [atom, atom2, atom3].index('CA')
            except:
                return 0
            cx_idx = [0, 1, 2]
            cx_idx.remove(n_idx)
            cx_idx.remove(ca_idx)
            cx_idx = cx_idx[0]
            if [atom, atom2, atom3][cx_idx][0] != 'C':
                return 0
            #print(f'{spec_type}: {w1} {w2} {w3} {seqidx} {seqidx2} {seqidx3} {atom} {atom2} {atom3} 1')
            return 1

        # N - CX - C    N(i)-CO(i-1)-CA/B(i-1) done, I dont have spectrum to try
        if spec_type == 'NCOCACB':
            try:
                c_idx = [w1, w2, w3].index(max([w1, w2, w3]))
            except:
                return 0
            try:
                cacb_idx = [w1, w2, w3].index(min([w1, w2, w3]))
            except:
                return 0
            n_idx = [0, 1, 2]
            n_idx.remove(c_idx)
            n_idx.remove(cacb_idx)
            n_idx = n_idx[0]
            n_seqidx = [seqidx, seqidx2, seqidx3][n_idx]
            c_seqidx = [seqidx, seqidx2, seqidx3][c_idx]
            cacb_seqidx = [seqidx, seqidx2, seqidx3][cacb_idx]
            if n_seqidx != c_seqidx + 1 or n_seqidx!=cacb_seqidx + 1:
                return 0
            if [atom, atom2, atom3][cacb_idx] not in ['CA', 'CB']:
                return 0
            if [atom, atom2, atom3][c_idx] != 'C':
                return 0
            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            return 1

        # N(i)-CO(i-1)-CX/C(i-1) /      Done
        if spec_type == 'NCOCX':
            try:
                n_idx = [atom, atom2, atom3].index('N')
            except:
                return 0
            try:
                c_idx = [w1, w2, w3].index(max([w1, w2, w3]))
            except:
                return 0
            cx_idx = [0, 1, 2]
            cx_idx.remove(n_idx)
            cx_idx.remove(c_idx)
            cx_idx = cx_idx[0]
            
            if [atom, atom2, atom3][cx_idx][0] != 'C':
                return 0
            cx_seqidx = [seqidx, seqidx2, seqidx3][cx_idx]
            n_seqidx = [seqidx, seqidx2, seqidx3][n_idx]
            c_seqidx = [seqidx, seqidx2, seqidx3][c_idx]
            
            if (n_seqidx - c_seqidx) != 1 or c_seqidx != cx_seqidx:
                return 0          
            if [atom, atom2, atom3][cx_idx] not in ['CA', 'CB', 'CG', 'CG1','CG2', 'CD', 'CD1', 'CD2', 'CE', 'C']:
                return 0

            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            #print(f'{spec_type}: {w1} {w2} {w3} {seqidx} {seqidx2} {seqidx3} {atom} {atom2} {atom3} 1')
            return 1
        
        # CA(i)-N(i)-CX/O(i-1)   Done
        if spec_type == 'CAN(CO)CX':
            try:    
                ca_idx = [atom, atom2, atom3].index('CA')
            except:
                return 0
            try:
                n_idx = [atom, atom2, atom3].index('N')
            except:
                return 0
            cx_idx = [0, 1, 2]
            cx_idx.remove(n_idx)
            cx_idx.remove(ca_idx)
            cx_idx = cx_idx[0]
            if [atom, atom2, atom3][cx_idx][0] != 'C':
                return 0
            n_seqidx = [seqidx, seqidx2, seqidx3][n_idx]
            ca_seqidx = [seqidx, seqidx2, seqidx3][ca_idx]
            cx_seqidx = [seqidx, seqidx2, seqidx3][cx_idx]
            
            if n_seqidx != ca_seqidx or (n_seqidx - cx_seqidx) != 1:
                return 0
            if [atom, atom2, atom3][ca_idx] != 'CA':
                return 0
            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            return 1
        
        # N-CA-C   N(i)-CA(i)-CO(i)  Done Not spectrum to try
        if spec_type == 'NCACO':
            if seqidx != seqidx3 or seqidx != seqidx2:
                return 0
            try:
                c_idx = [w1, w2, w3].index(max([w1, w2, w3]))
            except:
                return 0
            try:
                ca_idx = [w1, w2, w3].index(min([w1, w2, w3]))
            except:
                return 0
            n_idx = [0, 1, 2]
            n_idx.remove(c_idx)
            n_idx.remove(ca_idx)
            n_idx = n_idx[0]

            if [atom, atom2, atom3][ca_idx] != ['CA']:
                return 0
            if [atom, atom2, atom3][c_idx] != 'C':
                return 0
            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            return 1
        
        #CA(i) - N(i) - CO(i-1)  Done
        if spec_type == 'CANCO':
            try:
                c_idx = [w1, w2, w3].index(max([w1, w2, w3]))
            except:
                return 0
            try:
                ca_idx = [w1, w2, w3].index(min([w1, w2, w3]))
            except:
                return 0
            n_idx = [0, 1, 2]
            n_idx.remove(c_idx)
            n_idx.remove(ca_idx)
            n_idx = n_idx[0]
            
            n_seqidx = [seqidx, seqidx2, seqidx3][n_idx]
            ca_seqidx = [seqidx, seqidx2, seqidx3][ca_idx]
            c_seqidx = [seqidx, seqidx2, seqidx3][c_idx]
            
            if n_seqidx != ca_seqidx or (n_seqidx-c_seqidx) != 1 :
                return 0
            if [atom, atom2, atom3][ca_idx] != ['CA']:
                return 0
            if [atom, atom2, atom3][c_idx] != 'C':
                return 0
            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            return 1
        
        #N(i)-CO(i-1)-CA(i-1) Done
        if spec_type == 'NCOCA':
            try:
                c_idx = [w1, w2, w3].index(max([w1, w2, w3]))
            except:
                return 0
            try:
                ca_idx = [w1, w2, w3].index(min([w1, w2, w3]))
            except:
                return 0
            n_idx = [0, 1, 2]
            n_idx.remove(c_idx)
            n_idx.remove(ca_idx)
            n_idx = n_idx[0]
            
            n_seqidx = [seqidx, seqidx2, seqidx3][n_idx]
            ca_seqidx = [seqidx, seqidx2, seqidx3][ca_idx]
            c_seqidx = [seqidx, seqidx2, seqidx3][c_idx]
            
            if (n_seqidx-ca_seqidx) !=1 or (n_seqidx-c_seqidx) != 1 :
                return 0
            if [atom, atom2, atom3][ca_idx] != 'CA':
                return 0
            if [atom, atom2, atom3][c_idx] != 'C':
                return 0
            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            return 1
        
        #  CA(i)-N(i)-CA/O(i-1)  Done
        if spec_type == 'CAN(CO)CA':
            try:    
                n_idx = [atom, atom2, atom3].index('N')
            except:
                return 0
            try:
                ca_idx = [atom, atom2, atom3].index('CA')
            except:
                return 0
            cx_idx = [0, 1, 2]
            cx_idx.remove(n_idx)
            cx_idx.remove(ca_idx)
            cx_idx = cx_idx[0]
            if [atom, atom2, atom3][cx_idx] not in ['C','CA']:
                return 0
            n_seqidx = [seqidx, seqidx2, seqidx3][n_idx]
            ca_seqidx = [seqidx, seqidx2, seqidx3][ca_idx]
            cx_seqidx = [seqidx, seqidx2, seqidx3][cx_idx] 
            
            if n_seqidx != ca_seqidx or (n_seqidx-cx_seqidx) !=1:
                return 0
            if [atom, atom2, atom3][ca_idx] != 'CA':
                return 0
            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            return 1
        
        #CA(i)-N(i)-CO/A/B(i-1)  Done
        if spec_type == 'CAN(CO)CACB':
            try:    
                n_idx = [atom, atom2, atom3].index('N')
            except:
                return 0
            try:
                ca_idx = [atom, atom2, atom3].index('CA')
            except:
                return 0
            cx_idx = [0, 1, 2]
            cx_idx.remove(n_idx)
            cx_idx.remove(ca_idx)
            cx_idx = cx_idx[0]
            if [atom, atom2, atom3][cx_idx] not in ['C','CA','CB']:
                return 0
            n_seqidx = [seqidx, seqidx2, seqidx3][n_idx]
            ca_seqidx = [seqidx, seqidx2, seqidx3][ca_idx]
            cx_seqidx = [seqidx, seqidx2, seqidx3][cx_idx] 
            
            if n_seqidx != ca_seqidx or (n_seqidx-cx_seqidx) !=1:
                return 0
            if [atom, atom2, atom3][ca_idx] != 'CA':
                return 0
            if [atom, atom2, atom3][n_idx] != 'N':
                return 0
            return 1

    # -----------------------------------------------------------------------------
    #
    def add_resonance_to_list(self, shift, resonance, seq, seqidx, atom, prob):  
        if prob < 0.001:
            return
        if shift > 0.001 and shift < 200.0:
            if shift < 14.0:
                iShiftIdx = int((shift + 5.05) * 10.0)  # proton rule
            else:
                iShiftIdx = int(shift + 200.5)

            self.resonance_list[iShiftIdx -
                                1].append([shift, resonance, seq, seqidx, atom, prob])
            self.resonance_list[iShiftIdx].append(
                [shift, resonance, seq, seqidx, atom, prob])
            self.resonance_list[iShiftIdx +  
                                1].append([shift, resonance, seq, seqidx, atom, prob])
    # -----------------------------------------------------------------------------
    #

    def sort_resonance_list(self):
        for i in range(len(self.resonance_list)):
            temp_list = self.resonance_list[i]
            if len(temp_list) != 0:
                self.resonance_ascend_list[i] = sorted(
                    temp_list, key=getFirstKey)
                self.resonance_descend_list[i] = sorted(
                    temp_list, key=getFirstKey, reverse=True)

    # -----------------------------------------------------------------------------
    #
    def find_probables(self, peak, spec_type):
        probable_list = []

        res_list = [[], [], [], []]

        tol_list = []
        
        # fill res_list
        for dim in range(len(peak.frequency)):
            shift = peak.frequency[dim]

            if shift < 14.0:
                iShiftIdx = int((shift + 5.05) * 10.0)  # proton rule
            else:
                iShiftIdx = int(shift + 200.5)

            # get tolerances
            try:
                # if shift < 14.0:
                #     tol = float(self.htol_ppm.get())  #you get the tol from the user
                if shift < 80.0:
                    tol = float(self.ctol_ppm.get())
                elif shift < 140.0:
                    tol = float(self.ntol_ppm.get())
                else:
                    tol = float(self.ctol_ppm.get())
            except:
                if shift < 14.0:
                    tol = 0.02
                else:
                    tol = 0.2

            from_shift = shift-tol
            to_shift = shift+tol
            tol_list.append(tol)

            # iShiftIdx-1 -> go descending
            desc_list = self.resonance_descend_list[iShiftIdx-1]
            shift_list = self.resonance_ascend_list[iShiftIdx]
            asce_list = self.resonance_ascend_list[iShiftIdx+1]
            for [shift, resonance, seq, seqidx, atom, prob] in desc_list:
                if shift < from_shift:
                    break
                if shift >= from_shift and shift <= to_shift:
                    iExist = 0
                    for (shift2, resonance2, seq2, seqidx2, atom2, prob2) in res_list[dim]:
                        if shift == shift2 and resonance == resonance2 and seqidx == seqidx2 and atom == atom2:
                            iExist = 1
                            break
                    if iExist == 0:
                        res_list[dim].append(
                            [shift, resonance, seq, seqidx, atom, prob])
            # iShiftIdx -> go ascending
            for [shift, resonance, seq, seqidx, atom, prob] in shift_list:
                if shift > to_shift:
                    break
                if shift >= from_shift and shift <= to_shift:
                    iExist = 0
                    for (shift2, resonance2, seq2, seqidx2, atom2, prob2) in res_list[dim]:
                        if shift == shift2 and resonance == resonance2 and seqidx == seqidx2 and atom == atom2:
                            iExist = 1
                            break
                    if iExist == 0:
                        res_list[dim].append(
                            [shift, resonance, seq, seqidx, atom, prob])
            # iShiftIdx -> go ascending
            for [shift, resonance, seq, seqidx, atom, prob] in asce_list:
                if shift > to_shift:
                    break
                if shift >= from_shift and shift <= to_shift:
                    iExist = 0
                    for (shift2, resonance2, seq2, seqidx2, atom2, prob2) in res_list[dim]:
                        if shift == shift2 and resonance == resonance2 and seqidx == seqidx2 and atom == atom2:
                            iExist = 1
                            break
                    if iExist == 0:
                        res_list[dim].append(
                            [shift, resonance, seq, seqidx, atom, prob])

        # make probable_list, and check okay pattern
        for w1 in res_list[0]:
            for w2 in res_list[1]:
                sspine_diff = abs(w1[0]-peak.frequency[0])/tol_list[0] + \
                    abs(w2[0]-peak.frequency[1])/tol_list[1]
                if len(peak.frequency) == 2:  #2 dimensions
                    if self.check_pattern(spec_type, w1, w2, None, None) == 1:
                        atm_list = [w1[4], w2[4]]
                        seqidx_list = [w1[3], w2[3]]
                        prob_list = [w1[5], w2[5]]
                        try:
                            idx1 = atm_list.index('N')
                            idx2 = atm_list.index('CA')
                            if seqidx_list[idx1] == seqidx_list[idx2]:
                                # in this case, they need to be the same probability
                                if prob_list[idx1] != prob_list[idx2]:
                                    continue
                                combined_prob = w1[5]
                            else:
                                combined_prob = w1[5]*w2[5]
                        
                        except:
                            combined_prob = w1[5]*w2[5]
                        sspinelabel = '%s-%s[%.3f]:(%.3f,%.3f)' % (
                            w1[1], w2[1], combined_prob, peak.frequency[0], peak.frequency[1])
                        asgnlabel = '%s-%s[' % (w1[1], w2[1])
                        iFound = 0
                        for (prob, plabel, pdiff) in probable_list:
                            if plabel.find(asgnlabel) != -1:
                                iFound = 1
                                break
                        if iFound == 0:
                            probable_list.append(
                                [combined_prob, sspinelabel, sspine_diff])
                        continue

                if len(peak.frequency) == 3: #3 dimension
                    for w3 in res_list[2]:
                        sspine_diff = abs(w1[0]-peak.frequency[0])/tol_list[0] + abs(
                            w2[0]-peak.frequency[1])/tol_list[1] + abs(w3[0]-peak.frequency[2])/tol_list[2]
                        if self.check_pattern(spec_type, w1, w2, w3, None) == 1:
                            atm_list = [w1[4], w2[4], w3[4]]
                            seqidx_list = [w1[3], w2[3], w3[3]]
                            prob_list = [w1[5], w2[5], w3[5]]
                            try:
                                idx1 = atm_list.index('N')
                                idx2 = atm_list.index('CA')
                                if seqidx_list[idx1] == seqidx_list[idx2]:
                                    # in this case, they need to be the same probability
                                    if prob_list[idx1] != prob_list[idx2]:
                                        continue
                                    if prob_list[idx1] != 0:
                                        combined_prob = w1[5]*w2[5] * \
                                            w3[5] / prob_list[idx1]
                                    else:
                                        combined_prob = 0
                                else:
                                    combined_prob = w1[5]*w2[5]*w3[5]
                            except:
                                combined_prob = w1[5]*w2[5]*w3[5]

                            sspinelabel = '%s-%s-%s[%.3f]:(%.3f,%.3f,%.3f)' % (
                                w1[1], w2[1], w3[1], combined_prob, peak.frequency[0], peak.frequency[1], peak.frequency[2])
                            asgnlabel = '%s-%s-%s[' % (w1[1], w2[1], w3[1])
                            iFound = 0
                            for (prob, plabel, pdiff) in probable_list:
                                if plabel.find(asgnlabel) != -1:
                                    iFound = 1
                                    break
                            if iFound == 0:
                                probable_list.append(
                                    [combined_prob, sspinelabel, sspine_diff])
                            continue

                if len(peak.frequency) == 4:
                    for w3 in res_list[2]:
                        for w4 in res_list[3]:
                            sspine_diff = abs(w1[0]-peak.frequency[0])/tol_list[0] + abs(w2[0]-peak.frequency[1])/tol_list[1] + abs(
                                w3[0]-peak.frequency[2])/tol_list[2] + abs(w4[0]-peak.frequency[3])/tol_list[3]
                        if self.check_pattern(spec_type, w1, w2, w3, w4) == 1:
                            atm_list = [w1[4], w2[4], w3[4], w4[4]]
                            seqidx_list = [w1[3], w2[3], w3[3], w4[3]]
                            prob_list = [w1[5], w2[5], w3[5], w4[5]]
                            try:
                                idx1 = atm_list.index('N')
                                idx2 = atm_list.index('CA')
                                if seqidx_list[idx1] == seqidx_list[idx2]:
                                    # in this case, they need to be the same probability
                                    if prob_list[idx1] != prob_list[idx2]:
                                        continue
                                    if prob_list[idx1] != 0:
                                        combined_prob = w1[5]*w2[5] * \
                                            w3[5]*w4[5] / prob_list[idx1]
                                    else:
                                        combined_prob = 0
                                else:
                                    combined_prob = w1[5]*w2[5]*w3[5]*w4[5]
                            except:
                                combined_prob = w1[5]*w2[5]*w3[5]*w4[5]

                            sspinelabel = '%s-%s-%s-%s[%.3f]:(%.3f,%.3f,%.3f,%.3f)' % (
                                w1[1], w2[1], w3[1], w4[1], combined_prob, peak.frequency[0], peak.frequency[1], peak.frequency[2], peak.frequency[3])
                            asgnlabel = '%s-%s-%s-%s[' % (w1[1],
                                                          w2[1], w3[1], w4[1])
                            iFound = 0
                            for (prob, plabel, pdiff) in probable_list:
                                if plabel.find(asgnlabel) != -1:
                                    iFound = 1
                                    break
                            if iFound == 0:
                                probable_list.append(
                                    [combined_prob, sspinelabel, sspine_diff])
                            continue

        # python sort is not working for unknown reason. so we do bubble sort
        temp_list = sorted(probable_list, key=lambda probable: (
            probable[0], probable[2]), reverse=True)
        return temp_list


    # -----------------------------------------------------------------------------
    #

    # ssPINE-POKY module
    def create_sspine_labels(self, sspine_bb_lines, sspine_sc_lines):         
        # remove all floating ssPINE labels first
        for peak in self.session.selected_peaks():
            peak.selected = 0

        for (s_type, spec) in self.listbox.line_data: # this is the information added for the user spec type and spectrum
            for pLabel in spec.label_list():   
                if pLabel.peak == None:
                    if pLabel.text.find(']:(') != -1:
                        pLabel.selected = 1
                else:
                    pLabel.selected = 0
        self.session.command_characters("")

        # parse chemical shifts
        szSeq = ' '
        self.resonance_list = []  
        self.resonance_ascend_list = []
        self.resonance_descend_list = []
        for i in range(500):
            self.resonance_list.append([])
            self.resonance_ascend_list.append([])
            self.resonance_descend_list.append([])

        for line in sspine_bb_lines:
            try:
                splited = line.strip().split()
                iSeqIdx = int(splited[0])
                szA = myseq.aaa2a(splited[1])
                szSeq = szSeq + szA
                # read bb prob
                # Prob 1
                prob = float(splited[2])

                # N
                resonance = '%s%dN' % (szA, iSeqIdx)
                shift = float(splited[3])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'N', prob)

                # CO
                resonance = '%s%dC' % (szA, iSeqIdx)
                shift = float(splited[4])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'C', prob)

                # CA
                resonance = '%s%dCA' % (szA, iSeqIdx)
                shift = float(splited[5])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'CA', prob)

                # CB
                resonance = '%s%dCB' % (szA, iSeqIdx)
                shift = float(splited[6])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'CB', prob)

                # Prob 2
                prob = float(splited[7])

                # N
                resonance = '%s%dN' % (szA, iSeqIdx)
                shift = float(splited[8])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'N', prob)

                # CO
                resonance = '%s%dC' % (szA, iSeqIdx)
                shift = float(splited[9])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'C', prob)


                # Prob 3
                prob = float(splited[10])
                # N
                resonance = '%s%dN' % (szA, iSeqIdx)
                shift = float(splited[11])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'N', prob)

                # CO
                resonance = '%s%dC' % (szA, iSeqIdx)
                shift = float(splited[12])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'C', prob)

                # Prob 4
                prob = float(splited[13])
                # N
                resonance = '%s%dN' % (szA, iSeqIdx)
                shift = float(splited[14])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'N', prob)

                # CO
                resonance = '%s%dC' % (szA, iSeqIdx)
                shift = float(splited[15])
                self.add_resonance_to_list(
                    shift, resonance, szA, iSeqIdx, 'C', prob)
            except:
                continue

        # side chains
        for line in sspine_sc_lines:
            try:
                splited = line.strip().split()
                iSeqIdx = int(splited[1])
                szA = myseq.aaa2a(splited[2])
                szAtom = splited[3]
                resonance = '%s%d%s' % (szA, iSeqIdx, szAtom)
            except:
                self.session.show_message('ssPINE label error', splited)
                continue

            try:
                for i in range(0, 4):
                    prob = float(splited[4+2*i])
                    shift = float(splited[5+2*i])
                    self.add_resonance_to_list(
                        shift, resonance, szA, iSeqIdx, szAtom, prob)
            except:
                self.session.show_message('ssPINE label error2', splited)
                continue

        try:
            self.sort_resonance_list()   
        except:
            self.session.show_message('ssPINE label error', 'sort error')
            return 0

        # loop spectrum
        # generate labels and .pl files  
        sspinelabel_color = ['turquoise', 'orange', 'magenta', 'tomato', 'gray']
        for (s_type, spec) in self.listbox.line_data:
            try:
                f = open(spec.save_path+'_pl', 'w') 
            except:
                self.session.show_message(
                    'ssPINE label error', spec.save_path+'_pl is not writable.')
                return 0

            spec_probable_list = []
            for peak in spec.peak_list():
                probable_list = self.find_probables(peak, s_type)
                spec_probable_list.extend(probable_list)
                for i in range(len(probable_list)):
                    prob, sspinelabel, sspinediff = probable_list[i]
                    # make labels- turquoise, orange, margenta, tomato, gray-
                    if i < len(sspinelabel_color):
                        sspinecolor = sspinelabel_color[i]
                    else:
                        sspinecolor = 'gray'
                    poky.Label(spec, sspinelabel, sspinecolor, peak.frequency)
                    # write to file
                    # f.write(sspinelabel+'\n')
            sorted_probable_list = sorted(
                spec_probable_list, key=lambda prob: (prob[0], prob[2]), reverse=True)
            for prob, sspinelabel, sspinediff in sorted_probable_list:
                f.write(sspinelabel+'\n')  #here is writing the _pl file
            f.close()
        return 1
    # -----------------------------------------------------------------------------

    def assign_the_best(self):

        iTotal = 0
        iAssigned = 0
        iDuplicate = 0
        session = self.session

        # we will deselect all peaks first because we will delete floating labels at the end
        selected_peaks = session.selected_peaks()
        for peak in selected_peaks:
            peak.selected = 0

        # self.notice_text.set("* Check peak counts if PINE keeps failing (Use 'pv' tool).")
        for (s_type, spec) in self.listbox.line_data:
            if spec == None:
                continue
            if len(spec.selected_peaks()) == 0:
                peaks = spec.peak_list()
            else:
                peaks = spec.selected_peaks()

            for pLabel in spec.label_list():
                pLabel.selected = 0

            pParsedLabelList = []  # parse first. -10.0 (0) to 200.0 (2100).
            pLabelList = []
            pFloatList = []
            pProbList = []

            # generate empty list
            for i in range(2100):
                pParsedLabelList.append([])

            # parse labels
            for pLabel in spec.label_list():
                # M1CA-M1N-M1H[1.0]:(53.2120,124.321,8.23400)
                splited = pLabel.text.split(':')
                if len(splited) == 1:
                    # if it doesn't have position property, then continue.
                    continue

                szLabel = (splited[0].split('['))[0]
                szProb = (((splited[0].split('['))[1]).split(']'))[0]
                try:
                    fProb = float(szProb)
                    if fProb < 0.5:
                        continue        # probability lower than 0.5 will not be used
                except:
                    continue
                splited2 = splited[1].split(
                    '(')        # (53.2120,124.321,8.23400)
                splitedsplited = splited2[1].split(')')
                # 53.2120,124.321,8.23400
                pos = splitedsplited[0].split(',')

                if len(pos) == 0:
                    continue
                if pos[0] == '?':
                    continue
                try:
                    idx = int(eval(pos[0]) * 10.0 + 100.0)
                    if idx < 0 or idx > 2099:
                        continue
                    fPos = []
                    for p in pos:
                        fPos.append(eval(p))
                except:
                    continue
                pParsedLabelList[idx].append([fPos, fProb, szLabel, pLabel])

            iTotal = iTotal + len(peaks)

            for pPeak in peaks:
                if pPeak.is_assigned == 1:
                    continue
                pPeak.selected = 0
                del pLabelList[:]
                del pFloatList[:]
                del pProbList[:]

                idx = int(pPeak.frequency[0]*10.0 + 100.0)
                if idx < 1 or idx > 2098:
                    continue

                pTempLabelList = pParsedLabelList[idx-1] + \
                    pParsedLabelList[idx]+pParsedLabelList[idx+1]

                for (fPos, fProb, szLabel, pLabel) in pTempLabelList:
                    bSame = 1
                    for i in range(len(fPos)):
                        if abs(pPeak.frequency[i]*1000.0 - fPos[i]*1000.0) > 1:
                            bSame = 0
                            break
                    if bSame == 0:
                        continue
                    pFloatList.append(pLabel)
                    pLabelList.append(szLabel)
                    pProbList.append(fProb)

                # The best is the first
                if len(pFloatList) == 0:
                    continue

                # Get maximum probability
                pmax = 0.0
                idx = -1
                for i in range(len(pProbList)):
                    if pmax < pProbList[i]:
                        pmax = pProbList[i]
                        idx = i

                line_data = pLabelList[idx]
                pList = line_data.split('-')   # (M1CA, M1N, M1H)
                for i in range(len(pList)):
                    tAssign = sputil.split_group_atom(
                        pList[i])        # (M1, CA)
                    pPeak.assign(i, tAssign[0], tAssign[1])
                pPeak.show_assignment_label()

                for pLabel in pFloatList:
                    pLabel.selected = 1
                iAssigned = iAssigned + 1

            for i in range(len(peaks)):
                pks1 = peaks[i]
                if pks1.is_assigned != 1:
                    continue
                for j in range(len(peaks)):
                    if i == j:
                        continue
                    pks2 = peaks[j]
                    if pks2.is_assigned != 1:
                        continue
                    if pks1.assignment == pks2.assignment:
                        iDuplicate = iDuplicate + 1

        session.command_characters("")

        szAssigned = '%d of %d peaks assigned. %d peaks have duplicated assignments in the same spectrum.' % (
            iAssigned, iTotal, iDuplicate)
        self.session.show_message('Assign the Best by ssPINE', szAssigned)

    # -----------------------------------------------------------------------------
    #
    def web_cb(self):
        self.session.open_url('https://poky.clas.ucdenver.edu/ssPINE/')
    # -----------------------------------------------------------------------------
    #
    def user_group_cb(self):
        url = 'https://groups.google.com/g/pinenmr-user-group'
        import webbrowser
        webbrowser.open(url)

# -----------------------------------------------------------------------------
#
def run_sspine(session):
    d = sputil.the_dialog(sspine_dialog, session)
    d.sp.update_combobox_cb()
    d.show_window(1)
    return d
# -----------------------------------------------------------------------------
# Launched from POKY Notepad
if __name__.startswith('py_'):
  import __main__
  s = __main__.main_session  
  run_sspine(s)
## -----------------------------------------------------------------------------
##
##def run_sspine(session):
##    session.run_resource_module('sspinenmr.py')
