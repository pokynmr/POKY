#
# Simple NMR-STAR 3.1 Generation Script.
# by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)
#
# To run this script:
#   In Poky Notepad,
#     File -> Run Python Module

import __main__
s = __main__.main_session
proj = s.project

print('\n\n\n------------------------------------------------------')
print('POKY Simple NMR-STAR 3.1 Generation Script')
print('by Woonghee Lee, Ph.D. (woonghee.lee@ucdenver.edu)')
print('Department of Chemistry, University of Colorado Denver')
print('------------------------------------------------------')

# Select condition
selected_cond_names = s.show_conditionselectiondialog(
                    'Conditions to export', 1)
selected_cond_list = selected_cond_names.split('\t')
if len(selected_cond_list) == 0:
    raise SystemExit
print(selected_cond_list)
    
# Ask inclusion of unused resonances
unused_resonance = s.show_message_yes_no('Unused resonances', 
    'Do you want to include unused resonances?')

# Ask path
str_path = s.save_filedialog('Save as...', 'Any (*);; NMR-STAR 3.1 (*.str)', '')
if str_path == '':
    raise SystemExit

# -----------------------------------------------------------------------------
#
import myseq
def GetBiomoleculeType(session, cond_list=None):
  atgcu = 0
  no_atgcu = 0
  for condition in session.project.condition_list():
    if cond_list != None:
      if not condition.name in cond_list:
        continue
    for resn in condition.resonance_list():
      aaa = myseq.a2aaa(resn.group.symbol)
      if aaa == 'XXX': continue
      if aaa in ['ADE', 'GUA', 'CYT', 'URA', 'THY']: 
          atgcu = atgcu + 1
      else: no_atgcu = no_atgcu + 1
  if no_atgcu == 0 and atgcu > 0: return 'n'    
  return 'p'
# -----------------------------------------------------------------------------
#
def GenerateNmrStar(session, path, cond_list=None, unused_resonance=1):
  groups = []
  szSeqHeader = """save_protein
   _Entity.Sf_category                       entity
   _Entity.Sf_framecode                      protein
   _Entity.Entry_ID                          .
   _Entity.ID                                1
   _Entity.BMRB_code                         .
   _Entity.Name                              protein
   _Entity.Type                              polymer
   _Entity.Polymer_common_type               .
   _Entity.Polymer_type                      .
   _Entity.Polymer_type_details              .
   _Entity.Polymer_strand_ID                 .

   loop_
      _Entity_comp_index.ID
      _Entity_comp_index.Comp_ID

"""    
  szHeader = """save_chem_shift_list_1
   _Assigned_chem_shift_list.Sf_category                   assigned_chemical_shifts
   _Assigned_chem_shift_list.Sf_framecode                  chem_shift_list_1
   _Assigned_chem_shift_list.Entry_ID                      .
   _Assigned_chem_shift_list.ID                            1
   _Assigned_chem_shift_list.Chem_shift_1H_err             .
   _Assigned_chem_shift_list.Chem_shift_13C_err            .
   _Assigned_chem_shift_list.Chem_shift_15N_err            .
   _Assigned_chem_shift_list.Chem_shift_31P_err            .
   _Assigned_chem_shift_list.Chem_shift_2H_err             .
   _Assigned_chem_shift_list.Chem_shift_19F_err            .
   _Assigned_chem_shift_list.Error_derivation_method       .
   _Assigned_chem_shift_list.Details                       .
   _Assigned_chem_shift_list.Text_data_format              .
   _Assigned_chem_shift_list.Text_data                     .

   loop_  
      _Atom_chem_shift.ID
      _Atom_chem_shift.Comp_index_ID
      _Atom_chem_shift.Seq_ID
      _Atom_chem_shift.Comp_ID
      _Atom_chem_shift.Atom_ID
      _Atom_chem_shift.Atom_type
      _Atom_chem_shift.Atom_isotope_number
      _Atom_chem_shift.Val
      _Atom_chem_shift.Val_err
      _Atom_chem_shift.Ambiguity_code
      _Atom_chem_shift.Assigned_chem_shift_list_ID
      
"""
  szFooter = "\n   stop_\nsave_\n"
  import pinelayout

  # we don't know this is protein or nucleic acids.
  # we will determine using resonances.
  bt = GetBiomoleculeType(session, cond_list)

  # Seq ID, Seq ID, Seq, Atom, Atom type, Isotop, CS, CS_Dev, 2
  for condition in session.project.condition_list():
    if cond_list != None:
      if not condition.name in cond_list:
        print('Skipped condition: ' + condition.name)
        continue
    for resonance in condition.resonance_list():
      if resonance.peak_count == 0 and unused_resonance == 0:
        print('Skipped resonance: ' + resonance.group + \
              ' in condition ' + condition.name)
        continue
      if bt == 'p':
        aaa = myseq.a2aaa(resonance.group.symbol)
      else:
        aaa = myseq.na_a2aaa(resonance.group.symbol)
      if aaa == 'XXX': continue
      
      atom_type = resonance.atom.name[0]
      if atom_type in ['Q', 'M']: 
        atom_type = 'H'
      isotope = 1
      if atom_type == 'C': 
        isotope = 13
      elif atom_type == 'N':
        isotope = 15
      elif atom_type == 'P':
        isotope = 31
            
      group = [resonance.group.number, aaa, resonance.atom.name, atom_type, 
                isotope, resonance.frequency, resonance.deviation]
      groups.append(group)
  groups.sort()
  
  group_list2 = myseq.BuildSequence(session, 3)
  group_list2 = myseq.CheckSequence(group_list2) # sort and fill the gap
  
  if myseq.GapInSequence(group_list2) == 1: # there is a gap
    session.show_message('Gap', 
      'There is a gap between your assignment.\n' + \
      'You need to specify a sequence file.')
    group_list2 = myseq.LoadSeq(session, group_list2) # this fills a gap
    if group_list2 == None:
      session.show_message('Sequence', 
        'Sequence information is not satisfactory.')
      return 0
     
  f = open(path, 'w')
  # add header
  f.write(szSeqHeader)
  for group2 in group_list2:
    szAAA = myseq.a2aaa(group2[1])
    if bt == 'n': 
      szAAA = myseq.aaa2na_aaa(szAAA)
    szWrite = '%9d %s\n' % (group2[0], szAAA)
    f.write(szWrite)
  f.write('%s\n' % (szFooter))
  f.write(szHeader)
  iLine = 1      

  # THIS AMBIGUITY CODE NEEDS ATTENTION BY USER!            
  for group in groups:
    iAmb = 1
    if group[0] == None:
      continue
    if group[1] in pinelayout.pseudo_layout:
      if group[2] in pinelayout.pseudo_layout[group[1]]:
        iAmb = 2
    if group[1] in pinelayout.meta_layout:
      if group[2] in pinelayout.meta_layout[group[1]]:
        iAmb = 2                  
                  
    try:
      szWrite = "%6d %4d %4d %5s %5s %2s %2d %8.3f %6.2f %d  1\n" % \
                (iLine, group[0], group[0], group[1], 
                  group[2], group[3], group[4], group[5], 
                  group[6], iAmb)
  
      iLine = iLine+1
      f.write(szWrite)
    except:
      session.show_message('error', group)
      pass      

  f.write(szFooter)
  f.close()
  return 1

r = GenerateNmrStar(s, str_path, selected_cond_list, unused_resonance)

if r == 1:
  if s.show_message_yes_no('Export to NMRSTAR 3.1', 
         'Saved. Do you want to open in POKY Notepad?') == 1:      
    s.open_notepad(str_path)
