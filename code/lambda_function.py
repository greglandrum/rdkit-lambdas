# Copyright (C) 2018 greg Landrum and T5 Informatics GmbH
# All Rights Reserved
import json,os,requests
from urllib import parse
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import rdBase
Chem.WrapLogs()

def _get_swagger_help(base_url,method):
    resp = requests.get(f'{base_url}/apispec_1.json')
    if resp.status_code != 200:
        return None,'*Error*: could not retrieve service definition'
    o = resp.json()
    ps = o['paths']
    if not method in ps:
        return None,f'could not find method {method} in service definition'
    params = ps[method]['post']['parameters']
    exclude = ('smiles','mol')
    fields = [dict(title=p['name'],value=p['description'],short=False) for p in params if p['name'] not in exclude]
    return fields,'ok'

def _process_options(text, options):
    for entry in text:
        p = entry.find('=')
        if p== -1:
            raise ValueError("Bad options: %s"%str(entry))
        k = entry[:p]
        v = entry[p+1:]
        options[k] = v

def _nciresolver(nm):
    # As of Nov 2020 "verify=False" is required here
    resp = requests.get(
        'https://cactus.nci.nih.gov/chemical/structure/%s/smiles' % nm,verify=False)
    smiles = None
    errmsg = None
    if resp.status_code == 404:
        errmsg = f'Name "{nm}" could not be resolved'
    elif resp.status_code != 200:
        errmsg = 'problems accessing resolver service'
    else:
        smiles = resp.text
    return smiles, errmsg


def _chemblresolver(nm):
    resp = requests.get(
        'https://www.ebi.ac.uk/chembl/api/data/molecule/%s.json' % nm)
    smiles = None
    errmsg = None
    if resp.status_code == 404:
        errmsg = f'Name "{nm}" could not be resolved'
    elif resp.status_code != 200:
        errmsg = 'problems accessing resolver service'
    else:
        d = json.loads(resp.text)
        smiles = d['molecule_structures']['canonical_smiles']
    return smiles, errmsg


def _resolve(nm):
    if nm.find('CHEMBL') == 0:
        return _chemblresolver(nm)
    else:
        return _nciresolver(nm)


def _worker(command, text):
    if not text or text=='status':
        return "using RDKit version %s"%rdBase.rdkitVersion
    if text == 'help':
        obj,errmsg = _get_swagger_help('https://mol-renderer2-dev.t5ix.io',
        '/to_img/mol.png')
        if obj is None:
            return {'text':errmsg, 'response_type':'ephemeral'}
        return { 'text':f'Help for command {command}',
        'attachments':[{'fields':obj}]}


    text = text.split(' ')
    if command in ('/smiles', '/depict'):
        smiles = text[0]
    elif command in '/resolve':
        nm = parse.quote(text[0])
        smiles, errmsg = _resolve(nm)
        if errmsg is not None:
            return {'text': errmsg, 'response_type': 'in_channel'}
    else:
        raise ValueError("Bad Command: %s" % command)

    if command != '/depict':
        options = {'w': 400, 'h': 300}
    else:
        options = {'w': 200, 'h': 150}

    if len(text) > 1:
        _process_options(text[1:], options)

    mol = Chem.MolFromSmiles(smiles)
    fields = []
    if int(options.get('sanitize', 1)):
        if not mol:
            raise ValueError(
                "Molecule could not be processed.")
        canon_smiles = Chem.MolToSmiles(mol, True)
        if command != '/depict':
            fields.append({'title': 'Canonical Smiles',
                           'value': canon_smiles, 'short': False})
            fields.append(
                {'title': 'MolWt', 'value': '%.2f' % rdMolDescriptors._CalcMolWt(mol), 'short': True})
            calc = rdMolDescriptors.Properties(
                ['CrippenClogP', 'lipinskiHBA', 'lipinskiHBD', 'NumRotatableBonds', 'tpsa'])
            ps = calc.ComputeProperties(mol)
            fields.append(
                {'title': 'MolLogP', 'value': '%.2f' % ps[0], 'short': True})
            fields.append(
                {'title': 'LipinskiHBA', 'value': int(ps[1]), 'short': True})
            fields.append(
                {'title': 'LipinskiHBD', 'value': int(ps[2]), 'short': True})
            fields.append({'title': 'NumRotatableBonds',
                           'value': int(ps[3]), 'short': True})
            fields.append(
                {'title': 'TPSA', 'value': '%.2f' % ps[4], 'short': True})
    else:
        canon_smiles = smiles
    qsmi = parse.quote(canon_smiles)

    url_3d = 'http://3dmol.csb.pitt.edu/viewer.html?url=https://mol-renderer2-dev.t5ix.io/to_3d/mol.sdf?smiles=%s&type=sdf&style=stick' % qsmi
    img_url = 'https://mol-renderer2-dev.t5ix.io/to_img/mol.png?smiles=%s' % qsmi
    for k in options:
        img_url += '&%s=%s' % (k, parse.quote(str(options[k])))
    # 2. Return a JSON payload
    # See https://api.slack.com/docs/formatting and
    # https://api.slack.com/docs/attachments to send richly formatted messages
    if command == '/depict':
        au = ''
        title = ''
        title_link = ''
    else:
        au = 'RDKit Depictor'
        title = 'Molecule Info'
        title_link = url_3d
    return {
        # Uncomment the line below for the response to be visible to everyone
        # 'response_type': 'in_channel',
        #
        'text': '',  # 'text provided: %s'%text,
        'response_type': 'in_channel',
        'attachments': [
            {
                'fallback': 'Rendering of %s' % smiles,
                # 'pretext': '%s'%(canon_smiles),
                'pretext': '%s' % (canon_smiles),
                'author_name': au,
                # 'author_link':'http://www.t5informatics.com'
                'title_link': title_link,
                # 'author_icon': 'http://flickr.com/icons/bobby.jpg',
                'title': title,
                'fields': fields,
                # 'fields': [
                #     {
                #         'title': 'Priority',
                #         'value': 'High',
                #         'short': False
                #     }
                # ],
                'image_url': img_url
                # 'thumb_url': 'http://example.com/path/to/thumb.png'
            }
        ]
    }


def lambda_handler(event, context):
    # Parse the parameters you need
    token = event.get('token', None) 
    command = event.get('command', None)
    channel = event.get('channel_name', None)
    user_name = event.get('user_name', None)
    text = event.get('text', None)
    team_id = event.get('team_id', None)

    valid_tokens = os.environ.get('VALID_TOKENS').split(',')
    if token not in valid_tokens:
        raise ValueError("bad token: %s" % token)
    # print("ACCESS  |COMMAND: %(command)s|TEAM: %(team_id)s|CHANNEL: %(channel)s|USER_NAME: %(user_name)s|" %
    #       locals(), file=sys.stderr, flush=True)
    return _worker(command, text)
