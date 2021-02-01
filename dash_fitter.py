# -*- coding: utf-8 -*-
import os, sys
import json
import pickle
import numpy as np
from scipy.integrate import odeint

# Dash modules
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State, Event
app = dash.Dash()

#from TPD import TPD_data
import TPD
TPD.CACHING = True

# Model
R = 8.314 # gas constant, [J/(mol*K)]
NA = 6.022e23 # Avogrado's number, [mol^-1]
a0 = 3.61E-10 # Angstrom
nZ0 = np.array([4/(2*np.sqrt(3)), 4/(2*np.sqrt(6)), 4/(2*np.sqrt(6))]) * a0**(-2)/NA

def desorption(theta, t, r, T0, Hdes):
    dtheta_dt = -theta * 1e13 * np.exp( -(Hdes[0]*1000 + Hdes[1]*1000*theta) / (R*(T0 + r*t)))
    return dtheta_dt

def specific_desorption(theta, t, T, Hdes):
    dtheta_dt = -theta * 1e13 * np.exp( -(Hdes[0]*1000 + Hdes[1]*1000*theta) / (R*(T)))
    return dtheta_dt

def gauss(T, sigma):
    return 1./(sigma*(2*np.pi)**(1./2))*np.exp(-1./2*(T/sigma)**2)

# Defined constants
colors = {
    'background': '#EAEDED',
    'text': '#2471A3'
    }

###########################################################
###              APP LAYOUT                             ###
###########################################################

app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    # --- --- --- --- ---
    # HEADER SECTION
    # --- --- --- --- ---
    html.H1('TPD fitter', style={'textAlign': 'center', 'color': colors['text']}),
    html.Hr(),

    # --- --- --- --- ---
    # INPUTS
    # --- --- --- --- ---
    html.H2('Constants', style={'textAlign': 'left', 'color': colors['text']}),
    html.Div(id='const 1',
             style={'display': 'none'},
             children=json.dumps(
               {'heat ramp': 2,
                'start temperature': 77,
                'end temperature': 300,
                'c1 state': 0,
                'c1 Hdes': [0, 0],
                'c2 state': 0,
                'c2 Hdes': [0, 0],
                'c3 state': 0,
                'c3 Hdes': [0, 0],
                }),
            ),

    html.Div([
        html.Div([
            # General inputs
            html.Label('Heat ramp'),
            dcc.Input(id='heat ramp',
                type='number',
                value=2,
                ),
            html.Label('Start temperature'),
            dcc.Input(id='start temperature',
                type='number',
                value=77,
                ),
            html.Label('End temperature'),
            dcc.Input(id='end temperature',
                type='number',
                value=300,
                ),
            dcc.RadioItems(id='radio_custom_heating',
                options=[{'label': 'Heating from data', 'value': 1},
                         {'label': 'Regular heating', 'value': 0}],
                value=0,
                ),
            ], className="three columns"),

        html.Div([
            # parameter 1
            html.Label('Site 1'),
            dcc.RadioItems(id='c1 state',
                options=[{'label': 'On', 'value': 1},
                         {'label': 'Off', 'value': 0}],
                value=0,
                labelStyle={'display': 'inline-block'},
                ),
            dcc.Input(id='c1 Hdes1',
                type='number',
                value=0,
                ),
            dcc.Input(id='c1 Hdes2',
                type='number',
                value=0,
                ),
            html.Label('Ratio'),
            dcc.Input(id='c1 ratio',
                type='number',
                value=100,
                ),
            ], className="three columns"),

        html.Div([
            # parameter 2
            html.Label('Site 2'),
            dcc.RadioItems(id='c2 state',
                options=[{'label': 'On', 'value': 1},
                         {'label': 'Off', 'value': 0}],
                value=0,
                labelStyle={'display': 'inline-block'},
                ),
            dcc.Input(id='c2 Hdes1',
                type='number',
                value=0,
                ),
            dcc.Input(id='c2 Hdes2',
                type='number',
                value=0,
                ),
            html.Label('Ratio'),
            dcc.Input(id='c2 ratio',
                type='number',
                value=100,
                ),
            ], className="three columns"),

        html.Div([
            # parameter 3
            html.Label('Site 3'),
            dcc.RadioItems(id='c3 state',
                options=[{'label': 'On', 'value': 1},
                         {'label': 'Off', 'value': 0}],
                value=0,
                labelStyle={'display': 'inline-block'},
                ),
            dcc.Input(id='c3 Hdes1',
                type='number',
                value=0,
                ),
            dcc.Input(id='c3 Hdes2',
                type='number',
                value=0,
                ),
            html.Label('Ratio'),
            dcc.Input(id='c3 ratio',
                type='number',
                value=100,
                ),
            ], className="three columns"),
        ], className="row"),

    # --- --- --- --- ---
    # USER DATA
    # --- --- --- --- ---
    html.H2('Data', style={'textAlign': 'left', 'color': colors['text']}),
    html.Div([
        html.Div([
            #html.Label(' '),
            dcc.RadioItems(id='data on/off',
                options=[{'label': 'On', 'value': 1},
                         {'label': 'Off', 'value': 0}],
                value=0,
                ),
            html.Button('Load', id='load data'),
            ], className="three columns"),
        html.Div([
            html.Label('Timestamp (cinfdata)'),
            dcc.Input(id='data id',
                type='text',
                value='2020-06-29 01:02:13',#2018-08-28 22:15:37
                ),
            ], className="three columns"),
        html.Div([
            html.Label('Sigma'),
            dcc.Input(id='sigma',
                type='number',
                value=10,
                ),
            ], className="three columns"),
        html.Div([
            html.Label('Scale'),
            dcc.Input(id='data scale',
                type='text',
                value='1e9',
                ),
            ], className="three columns"),
        ], className="row"),
    html.Div(id='data options', children=json.dumps([False, [], []]), style={'display': 'none'}),
    html.Div(id='frozen heat ramp', children=json.dumps([False, [], []]), style={'display': 'none'}),
    html.Div(id='user data selection', style={'display': 'none'},
        children=json.dumps(None)),
    html.Div(id='user data options', children=[
            dcc.Checklist(id='user data index',
                labelStyle={'display': 'inline-block'},
                options=[],
                values=[],
                style={'display': 'none'}
                ),
            dcc.Checklist(id='user data label',
                labelStyle={'display': 'inline-block'},
                options=[],
                values=[],
                style={'display': 'none'}
                ),
        ]),
    html.Hr(),
    
    # --- --- --- --- ---
    # FIGURES
    # --- --- --- --- ---
    html.Button('Update', id='update graph'),
    html.Button('Save fig 2', id='save_figure_2'),
    html.Div(id='dummy data', children=json.dumps(([0,1],[], [])), style={'display': 'none'}),
    html.Div([
        html.Div(
            dcc.Graph(id='tpd 1'),
            className="three columns"),
        html.Div(
            dcc.Graph(id='tpd 2'),
            className="six columns"),
        html.Div(
            dcc.Graph(id='tpd 3'),
            className="three columns"),
        ], className="row"),

])


# CSS style
app.css.append_css({
    'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
})
###########################################################
###              CALLBACKS                              ###
###########################################################

# HEADER
#--------------
@app.callback(Output('const 1', 'children'),
    [Input('heat ramp', 'value'),
    Input('start temperature', 'value'),
    Input('end temperature', 'value'),
    Input('c1 state', 'value'),
    Input('c1 Hdes1', 'value'),
    Input('c1 Hdes2', 'value'),
    Input('c1 ratio', 'value'),
    Input('c2 state', 'value'),
    Input('c2 Hdes1', 'value'),
    Input('c2 Hdes2', 'value'),
    Input('c2 ratio', 'value'),
    Input('c3 state', 'value'),
    Input('c3 Hdes1', 'value'),
    Input('c3 Hdes2', 'value'),
    Input('c3 ratio', 'value'),
    Input('radio_custom_heating', 'value'),
], state=[State('const 1', 'children')])
def update_settings(ramp, T0, T1, c1s, c1H1, c1H2, c1r, c2s, c2H1, c2H2, c2r, c3s, c3H1, c3H2, c3r, heating_type, data_string):
    data = json.loads(data_string)

    data['heat ramp'] = ramp
    data['start temperature'] = T0
    data['end temperature'] = T1
    data['c1 state'] = c1s
    data['c1 Hdes'] = [c1H1, c1H2]
    data['c1 ratio'] = c1r
    data['c2 state'] = c2s
    data['c2 Hdes'] = [c2H1, c2H2]
    data['c2 ratio'] = c2r
    data['c3 state'] = c3s
    data['c3 Hdes'] = [c3H1, c3H2]
    data['c3 ratio'] = c3r
    data['custom heating'] = heating_type

    return json.dumps(data)

@app.callback(Output('dummy data', 'children'),
    state=[State('const 1', 'children'),
           State('data options', 'children'),
           State('data id','value'),
],
    events=[Event('update graph', 'click')])
def compute_data(data_string, data, timestamp):
    options = json.loads(data_string)
    is_on, mass_labels, index = json.loads(data)
    data = TPD.Experiment(timestamp)

    if options['custom heating']:
        function = specific_desorption
        interpolate = interp1d(times, temps, kind='linear', bounds_error=False, fill_value='extrapolate')
        pass
    else:
        function = desorption
        t_total = (options['end temperature'] - options['start temperature'])/options['heat ramp']
        t_vector = np.linspace(0, t_total, t_total*100 + 1)
    T_vector = t_vector*options['heat ramp'] + options['start temperature']
    #t_vector = (options['end temperature'] - options['start temperature'])/options['heat ramp']

    if True: # Pause heating at -100 C = 173 K
        mask = np.where(T_vector > 173)[0]
        print(T_vector[mask[0]:mask[0]+10])
        T_vector[mask] = T_vector[mask[0]]
    """
def desorption(theta, t, r, T0, Hdes):
    dtheta_dt = -theta * 1e13 * np.exp( -(Hdes[0]*1000 + Hdes[1]*1000*theta) / (R*(T0 + r*t)))
    return dtheta_dt

def specific_desorption(theta, t, T, Hdes):
    dtheta_dt = -theta * 1e13 * np.exp( -(Hdes[0]*1000 + Hdes[1]*1000*theta) / (R*(T)))
    return dtheta_dt
    """

    theta = []
    if False:
        if options['c1 state']:
            theta1 = odeint(desorption, 1, t_vector, args=(options['heat ramp'], options['start temperature'], options['c1 Hdes']))
            theta1 = theta1.transpose()[0]
            theta.append(theta1.tolist())

        if options['c2 state']:
            theta2 = odeint(desorption, 1, t_vector, args=(options['heat ramp'], options['start temperature'], options['c2 Hdes']))
            theta2 = theta2.transpose()[0]
            theta.append(theta2.tolist())

        if options['c3 state']:
            theta3 = odeint(desorption, 1, t_vector, args=(options['heat ramp'], options['start temperature'], options['c3 Hdes']))
            theta3 = theta3.transpose()[0]
            theta.append(theta3.tolist())
    else:
        if options['c1 state']:
            theta1 = odeint(specific_desorption, np.ones(len(t_vector)), t_vector, args=(T_vector, options['c1 Hdes']))
            theta1 = theta1.transpose()[0]
            theta.append(theta1.tolist())

        if options['c2 state']:
            theta2 = odeint(specific_desorption, np.ones(len(t_vector)), t_vector, args=(T_vector, options['c2 Hdes']))
            theta2 = theta2.transpose()[0]
            theta.append(theta2.tolist())

        if options['c3 state']:
            theta3 = odeint(specific_desorption, np.ones(len(t_vector)), t_vector, args=(T_vector, options['c3 Hdes']))
            theta3 = theta3.transpose()[0]
            theta.append(theta3.tolist())

    return json.dumps([t_vector.tolist(), theta, T_vector.tolist()])

@app.callback(Output('data options', 'children'),
    state=[
        State('data id','value'),
        State('data on/off','value')],
    events=[
        Event('load data', 'click')
    ])
def load_data(timestamp, is_on):
    """Load data in cinfdata/omicron/timestamp into hidden Div element.
    """
    output = dict()
    #experiment = TPD_data(timestamp, caching=False)
    #experiments = experiment.isolate_experiments()
    experiment = TPD.Experiment(timestamp)
    experiments = experiment.exps
    #DATA = experiment

    # Pick out mass data
    mass_labels = []
    for label in experiment.labels:
        if label.startswith('M') and label[1] in [str(i) for i in range(10)]:
            mass_labels.append(label)
    mass_labels.sort()
    #print(dir(experiment))

    # Append data to ouput dict
    #output = dict()
    #print(mass_labels)
    index = list(experiments.keys())
    #print(index)
    #for label in mass_labels:
    #    output[label] = dict()
    #    for i in index:
    #        data = experiments[i][label]
    #        output[label][i] = [list(data[0]), list(data[1]-min(data[1])), list(data[2]+273.15), list(data[3])]
    return json.dumps([is_on, mass_labels, index])

@app.callback(
    Output('user data options','children'),
    [Input('data options', 'children')]
    )
def generate_udic(options):
    """Generate checklist of user data index
    """
    print('udic')
    options = json.loads(options)
    print(repr(options))
    data_on, labels, index = options
    children = []
    print(index)
    
    # Checklist 1: index
    cl1 = dcc.Checklist(id='user data index',
        labelStyle={'display': 'inline-block'},
        options=[{'label': key, 'value': key} for key in index],
        values=[],
        )
    children.append(cl1)

    # Checklist 2: labels
    cl2 = dcc.Checklist(id='user data label',
        labelStyle={'display': 'inline-block'},
        options=[{'label': label, 'value': label} for label in labels],
        values=[],
        )
    children.append(cl2)

    # Return elements
    return children

@app.callback(
    Output('user data selection', 'children'),
    [Input('user data index', 'values'),
    Input('user data label', 'values')]
)
def update_selection(index, label):
    """Send new information about user choice of index and label
    """
    print(index, label)
    return_string = json.dumps({'index': index, 'label': label})
    print(return_string)
    return return_string

# FIGURES
#--------------
@app.callback(Output('tpd 2', 'figure'),
    [Input('dummy data', 'children')],
    state=[State('const 1', 'children'),
           State('sigma', 'value'),
           State('data options', 'children'),
           State('user data index', 'values'),
           State('user data label', 'values'),
           State('data id','value'),

           State('data scale', 'value')])
def fig2(data_string, meta_string, sigma, user_data, uindex, ulabel, timestamp, scale):

    time, theta, temperature = json.loads(data_string)
    time = np.array(time)
    meta = json.loads(meta_string)
    data_on, labels, index = json.loads(user_data)
    data = TPD.Experiment(timestamp)
    print('Fig 2:')
    print(uindex)
    print(ulabel)
    print('---')
    ramp = meta['heat ramp']
    T0 = meta['start temperature']

    x = time*ramp + T0
    dT = x[1] - x[0]
    #sigma = 10
    T_vec = np.arange(-3*sigma, 3*sigma+1, dT)
    Gauss_function = gauss(T_vec, sigma)
    
    # weight
    Z = [meta['c{} ratio'.format(i+1)] for i in range(3)]

    # Data sets
    data_traces = []
    maximum = 0
    if data_on:
        print(data.exps.keys())
        # Append selected data to traces
        for mass in ulabel:
            #print(mass)
            #print(data.exps[i][mass].keys())
            for i in uindex:
                #temp = data[mass][str(i)]
                dat = data.exps[int(i)][mass]
                peak = max(dat[1])
                if peak > maximum:
                    maximum = peak
                data_traces.append(
           {#'x': dat[2]+273.15,
            'x': time,
            'y': dat[1],
            'type': 'line',
            'name': '{}-{}'.format(i, mass)
            })

    # Simulated sets
    traces = [
            {#'x': x,
             'x': time,
             'y': np.convolve(
                desorption(np.array(cov), time, ramp, T0, meta['c{} Hdes'.format(i+1)])*(-Z[i])*maximum,
                Gauss_function, mode='same'),
             'type': 'line',
             'name': 'Site {}'.format(i+1)}
            for i, cov in enumerate(theta)]
    sum_trace = np.zeros(len(x))
    for i, trace in enumerate(traces):
        sum_trace += trace['y']
    traces.append(
            {#'x': x,
             'x': time,
             'y': sum_trace,
             'type': 'line',
             'name': 'Total',
            })
    traces.append(data_traces)


    return {
        'data': traces,
        'layout': {
            'plot_bgcolor': colors['background'],
            'paper_bgcolor': colors['background'],
            'font': {'color': colors['text']},
            'title': 'Chosen data slope',
            'xaxis': dict(
                #title='Temperature (K)'),
                title='Time (s)'),
            'yaxis': dict(
                title='dT/dt'),
            },
        }

@app.callback(Output('tpd 1', 'figure'),
    [Input('dummy data', 'children')],
    state=[State('const 1', 'children'),
        State('data options', 'children'),
        State('data id','value'),
])
def fig1(data_string, meta_string, data_options, timestamp):

    time, theta, temperature = json.loads(data_string)
    time, theta, temperature = np.array(time), np.array(theta), np.array(temperature)
    print(time[0], time[-1])
    print(temperature[0], temperature[-1])
    print(theta[0], theta[-1])
    meta = json.loads(meta_string)
    data_on, labels, index = json.loads(data_options)
    data = TPD.Experiment(timestamp)

    ramp = meta['heat ramp']
    T0 = meta['start temperature']
    

    traces = [{#'x': time*ramp + T0,
             #'y': cov,
             'x': time,
             'y': temperature,
             'type': 'line',
             'name': 'Site {}'.format(i+1)}
            for i, cov in enumerate(theta)]


    return {
        'data': traces,
        'layout': {
            'plot_bgcolor': colors['background'],
            'paper_bgcolor': colors['background'],
            'font': {'color': colors['text']},
            'title': 'Site coverage of temperature',
            'xaxis': dict(
                #title='Temperature (K)'),
                title='Time (s)'),
            'yaxis': dict(
                #title='Coverage'),
                title='Temperature (K)'),
            },
        }

@app.callback(Output('tpd 3', 'figure'),
    [Input('dummy data', 'children')],
    state=[State('const 1', 'children'),
        State('data options', 'children'),
        State('data id','value'),
])
def fig3(data_string, meta_string, data_options, timestamp):

    time, theta, temperature = json.loads(data_string)
    time, theta, temperature = np.array(time), np.array(theta), np.array(temperature)
    print(time[0], time[-1])
    print(temperature[0], temperature[-1])
    print(theta[0], theta[-1])
    meta = json.loads(meta_string)
    data_on, labels, index = json.loads(data_options)
    #data = TPD.Experiment(timestamp)

    ramp = meta['heat ramp']
    T0 = meta['start temperature']

    traces = [{
             'x': time,
             'y': cov,
             'type': 'line',
             'name': 'Site {}'.format(i+1)}
            for i, cov in enumerate(theta)]


    return {
        'data': traces,
        'layout': {
            'plot_bgcolor': colors['background'],
            'paper_bgcolor': colors['background'],
            'font': {'color': colors['text']},
            'title': 'Site coverage of temperature',
            'xaxis': dict(
                #title='Temperature (K)'),
                title='Time (s)'),
            'yaxis': dict(
                title='Coverage'),
            },
        }


###########################################################
if __name__ == '__main__':
    app.run_server(debug=True)


