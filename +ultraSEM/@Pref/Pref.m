classdef Pref
%ULTRASEM.PREF   Class for managing preferences in the ULTRASEM system.
%   ULTRASEM.PREF is a class for managing ULTRASEM construction-time and
%   solver preferences.
%
%   Available preferences:
%
%     discretization      - Discretization type.
%       ['coeffs']
%        'values'
%
%       This option determines whether PDEs are discretized using the
%       ultraspherical spectral method ('coeffs') or spectral collocation
%       ('values'). (Not currently supported.)
%
%     discSize            - Default discretization size.
%       [21]
%
%        This option specifies the default discretization size to use on a
%        leaf when no discretization size is given by the user.
%
%     interfaceDegree     - Routine to reconcile different p at an interface.
%       [@max]
%        @min
%        @mean
%        @(p1,p2) ...
%
%        This option determines the polynomial degree that is used along an
%        interface when the elements on either side of the interface use
%        different polynomial degrees, p1 and p2. Any function handle of the
%        form @(p1,p2) ... that returns an integer can be specified.
%
%     solver              - Linear solver to use on leaves.
%        '\'
%       ['woodbury']
%        'LU'
%
%        This option specifies the linear solver to use when constructing the
%        solution operator on each leaf.
%
%     splitTriangles      - Should we split triangles into quadrilaterals?
%       [true]
%        false
%
%        If true, ULTRASEM.TRIANGLEs are split into three quadrilaterals. If
%        false, ULTRASEM.TRIANGLEs are mapped to [-1,1]^2 via the Duffy
%        transformation and then discretized.
%
%   Constructor inputs:
%
%       P = ULTRASEM.PREF() creates an ULTRASEM.PREF object with the default
%       values of the preferences. For a list of all available preferences,
%       see above.
%
%       P = ULTRASEM.PREF(Q), where Q is a MATLAB structure, uses the
%       field/value pairs in Q to set the properties of P. If a field of Q
%       has a name which matches a property of P, the value of that property
%       of P is set to the value associated to that field in Q. If a field of
%       Q has a name that does not correspond to a known preference, then an
%       error is thrown.
%
%       P = ULTRASEM.PREF(Q), where Q is an ULTRASEM.PREF, sets P to be a
%       copy of Q.

%   Copyright 2020 Dan Fortunato, Nick Hale, and Alex Townsend.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties ( Access = protected )

        % MATLAB struct to hold a list of preferences.
        prefList

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods

        function outPref = Pref(varargin)

            if ( nargin == 0 )
                Q = struct();
            elseif ( nargin == 1 )
                if ( isa(varargin{1}, 'ultraSEM.Pref') )
                    outPref = varargin{1};
                    return
                else
                    Q = varargin{1};
                end
            else
                error('ULTRASEM:PREF:Pref:inputs', ...
                    'Too many input arguments.');
            end

            % Initialize default preference values.
            outPref.prefList = ultraSEM.Pref.manageDefaultPrefs('get');

            % Copy fields from Q, merging incomplete substructures.
            for field = fieldnames(Q).'
                field1 = field{1};
                if ( isfield(outPref.prefList, field1) )
                    if ( isstruct(outPref.prefList.(field1)) )
                        outPref.prefList.(field1) = mergePrefStructs( ...
                            outPref.prefList.(field1),                ...
                            Q.(field1)                                ...
                        );
                    else
                        if ( ultraSEM.Pref.isValidPrefVal(field1, Q.(field1)) )
                            outPref.prefList.(field1) = Q.(field1);
                        else
                            error('ULTRASEM:PREF:Pref:invalidPrefVal', ...
                                'Invalid preference value.');
                        end
                    end
                else
                    error('ULTRASEM:PREF:Pref:badPref', ...
                        'Unrecognized preference name.');
                end
            end

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static )

        % Merge preference structures.
        pref = mergePrefStructs(pref1, pref2);

        % Set default preferences.
        setDefaults(varargin);

        % Get factory default preferences.
        pref = getFactoryDefaults();

        % Get valid preference values for a specific preference.
        prefVals = getValidPrefVals(prefName);

        % Determine if a preference setting is valid.
        isValid = isValidPrefVal(prefName, prefVal);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRIVATE STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    methods ( Static, Access = private )

        % Private method for handling default preferences.
        varargout = manageDefaultPrefs(varargin);

        % Get structure of factory default preferences.
        pref = factoryDefaultPrefs();

    end

end
