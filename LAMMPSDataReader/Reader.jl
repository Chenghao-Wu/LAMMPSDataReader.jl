export  Universe,
        box,
        mass,
        types,
        moleid,
        charges,
        positions,
        bonds,
        bondtype,
        angles,
        angletype,
        velocities

# Sections will all start with one of these words
# and run until the next section title
SECTIONS = [
    "Atoms",  # Molecular topology sections
    "Velocities",
    "Masses",
    "Ellipsoids",
    "Lines",
    "Triangles",
    "Bodies",
    "Bonds",  # Forcefield sections
    "Angles",
    "Dihedrals",
    "Impropers",
    "Pair",
    "Pair LJCoeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "BondBond Coeffs",  # Class 2 FF sections
    "BondAngle Coeffs",
    "MiddleBondTorsion Coeffs",
    "EndBondTorsion Coeffs",
    "AngleTorsion Coeffs",
    "AngleAngleTorsion Coeffs",
    "BondBond13 Coeffs",
    "AngleAngle Coeffs",
]

HEADERS = [
    "atoms",
    "bonds",
    "angles",
    "dihedrals",
    "impropers",
    "atom types",
    "bond types",
    "angle types",
    "dihedral types",
    "improper types",
    "extra bond per atom",
    "extra angle per atom",
    "extra dihedral per atom",
    "extra improper per atom",
    "extra special per atom",
    "ellipsoids",
    "lines",
    "triangles",
    "bodies",
    "xlo xhi",
    "ylo yhi",
    "zlo zhi",
    "xy xz yz",
]

function get_StartLineNumber(lines)
    starts=[]
    Sects=Dict{String,Any}()
    for (line_index,line) in enumerate(lines)
        line=split(line,"#")
        if length(line)>0
            #println(line[1])
            if strip(line[1]) in SECTIONS
                push!(starts,line_index)
            end
        end
    end
    for (i,sectii) in enumerate(starts[1:length(starts)-1])
        push!(Sects,strip(split(lines[sectii],"#")[1])=>lines[sectii+2:starts[i+1]])
    end
    Sects
end

function get_TokenLineNumber(lines)
    headers=Dict{String,Any}()
    for line in lines
        for token in HEADERS
            if endswith(line,token)
                push!(headers,token=>strip(split(line,token)[1]))
            end
        end
    end
    headers
end

function _parse_mass(dataline)
    masses = Dict{Int,Float64}()
    for line in dataline
        line = split(line)
        if length(line)==2
            push!(masses,parse(Int64, line[1])=>parse(Float64, line[2]))
        end
    end

    return masses
end

function _determine_linenumber(datalines)
    linenumber=0
    for line in datalines
        line = split(line)
        if length(line)>3
            linenumber+=1
        end
    end
    linenumber
end

function _parse_pos(datalines,style_dict)

    if style_dict == false
        if length(split(datalines[1])) in (7, 10)
            style_dict = Dict("id"=> 1, "x"=> 5, "y"=> 6, "z"=> 7)
        else
            style_dict =  Dict("id"=> 1, "x"=> 4, "y"=> 5, "z"=> 6)
        end
    else
        style_dict = style_dict
    end

    number_atoms=_determine_linenumber(datalines)

    poses=zeros(3,number_atoms)
    ids=zeros(Int64,number_atoms)

    for lineii=1:number_atoms
        line = split(datalines[lineii])
        if length(line)>3
            #println(style_dict["id"])
            ids[lineii]=parse(Int64,line[style_dict["id"]])
            poses[1,lineii]=parse(Float64, line[style_dict["x"]])
            poses[2,lineii]=parse(Float64, line[style_dict["y"]])
            poses[3,lineii]=parse(Float64, line[style_dict["z"]])
        end
    end

    order = sortperm(ids)
    poses = poses[order]

    return poses,order
end

function _parse_vel(datalines,order)

    number_atoms=_determine_linenumber(datalines)

    vels=zeros(3,number_atoms)
    
    for lineii=1:number_atoms
        line = split(datalines[lineii])
        vels[1,lineii]=parse(Float64, line[2])
        vels[2,lineii]=parse(Float64, line[3])
        vels[3,lineii]=parse(Float64, line[4])
        
    end

    vels = vels[order]

    return vels
end

function _parse_bond_section(datalines,nentries)
    """Read lines and strip information
    Arguments
    ---------
    datalines : list
        the raw lines from the data file
    nentries : int (2,3,4)=>(bond,angle,dihedral)
        number of integers per line
    mapping : dict
        converts atom_ids to index within topology
    Returns
    -------
    types : tuple of strings
        type of the bond/angle/dihedral,
    indices : tuple of ints
        indices of atoms involved
    """
    number_bonds=_determine_linenumber(datalines)

    section=[]
    type=[]
    for lineii=1:number_bonds
        line=split(datalines[lineii])
        push!(section,[parse(Int64,x) for x in line[3:3+nentries]])
        push!(type,parse(Int64,line[2]))
    end
    type,section
end

function _parse_atoms(datalines,style_dict,massdic)
    """Creates a Topology object
    Adds the following attributes
        - resid
        - type
        - masses (optional)
        - charge (optional)
    Lammps atoms can have lots of different formats,
    and even custom formats
    http://lammps.sandia.gov/doc/atom_style.html
    Treated here are
    - atoms with 7 fields (with charge) "full"
    - atoms with 6 fields (no charge) "molecular"
    Arguments
    ---------
    datalines - the relevent lines from the data file
    massdict - dictionary relating type to mass
    Returns
    -------
    top - Topology object
    """

    number_atoms=_determine_linenumber(datalines)

    sd=0

    if style_dict == false
        sd = Dict("id"=> 1, "resid"=> 2, "type"=> 3)
        n=length(split(datalines[1]))
        if n in (7,10)
            push!(sd , "charge"=> 4)
        end
    else
        sd = style_dict
    end

    has_charge = "charge" in keys(sd)
    has_resid = "resid" in keys(sd)

     # atom ids aren't necessarily sequential
     atom_ids = zeros(Int64,number_atoms)
     types = zeros(Int64,number_atoms)

     if has_resid
        resids = zeros(Int64,number_atoms)
    else
        resids = ones(Int64,number_atoms)
    end

    if has_charge
        charges = zeros(number_atoms)
    end

    for lineii=1:number_atoms
        line=split(datalines[lineii])
        # these numpy array are already typed correctly,
        # so just pass the raw strings
        # and let numpy handle the conversion
        atom_ids[lineii]=parse(Int64,line[sd["id"]])
        if has_resid
            resids[lineii]=parse(Int64,line[sd["resid"]])
        end
        types[lineii]=parse(Int64,line[sd["type"]])
        if has_charge
            charges[lineii]=parse(Float64,line[sd["charge"]])
        end
    end

    # at this point, we've read the atoms section,
    # but it's still (potentially) unordered
    # TODO: Maybe we can optimise by checking if we need to sort
    # ie `if np.any(np.diff(atom_ids) > 1)`  but we want to search
    # in a generatorish way, np.any() would check everything at once
    order = sortperm(atom_ids)
    atomids=atom_ids[order]
    
    types=types[order]

    masses=[massdic[typeii] for typeii in types]

    if has_resid
        resids=resids[order]
    end
    if has_charge
        charges=charges[order]
    end

    return types,atomids,charges,masses
end

function _parse_box(header)
    x1, x2 = (parse(Float64,x) for x in split(header["xlo xhi"]))
    x = x2 - x1
    y1, y2 = (parse(Float64,x) for x in split(header["ylo yhi"]))
    y = y2 - y1
    z1, z2 = (parse(Float64,x) for x in split(header["zlo zhi"]))
    z = z2 - z1

    unitcell = zeros(3)
    unitcell[1:3] .= x, y, z
    unitcell
end

struct Universe
    sects
    headers
end


function Universe(filename::String)
    file_in=open(filename)
    lines=readlines(file_in)
    sects=get_StartLineNumber(lines)
    headers=get_TokenLineNumber(lines)
    return Universe(sects,headers)
end

#try
#    masses=_parse_mass(sects["Masses"])
#catch e
#    @error "LAMMPSDataReader.jl can not find Masses for atoms"
#end
#
#if "Atoms" ∉ keys(sects)
#    @error "LAMMPSDataReader.jl can not find Atoms"
#end
#
#try 
#    atoms=_parse_atoms(sects["Atoms"],false,masses)
#catch e
#    @error "LAMMPSDataReader.jl Failed to parse atoms section."
#end

function box(universe)
    box=nothing
    try 
        box=_parse_box(universe.headers)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse box header."
    end
    box
end

function mass(universe)
    masses=_parse_mass(universe.sects["Masses"])
    atoms=nothing
    try 
        atoms=_parse_atoms(universe.sects["Atoms"],false,masses)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse mass section."
    end
    atoms[4]
end

function types(universe)
    masses=_parse_mass(universe.sects["Masses"])
    atoms=nothing
    try 
        atoms=_parse_atoms(universe.sects["Atoms"],false,masses)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse mass section."
    end
    atoms[1]
end

function moleid(universe)
    masses=_parse_mass(universe.sects["Masses"])
    atoms=nothing
    try 
        atoms=_parse_atoms(universe.sects["Atoms"],false,masses)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse mass section."
    end
    atoms[2]
end

function charges(universe)
    masses=_parse_mass(universe.sects["Masses"])
    atoms=nothing
    try 
        atoms=_parse_atoms(universe.sects["Atoms"],false,masses)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse mass section."
    end
    atoms[3]
end

function positions(universe)
    position=nothing
    try 
        position=_parse_pos(universe.sects["Atoms"],false)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse position section."
    end
    position[1]
end


function velocities(universe)
    velocity=nothing
    try 
        order=_parse_pos(universe.sects["Atoms"],false)[2]
        velocity=_parse_vel(universe.sects["Velocities"],order)
        
    catch KeyError
        @error "LAMMPSDataReader.jl Failed to parse velocities section."
    end
    velocity
end

function bonds(universe)
    bond=nothing
    try 
        bond=_parse_bond_section(universe.sects["Bonds"],1)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse bonds section."
    end
    bond[2]
end

function bondtype(universe)
    bond=nothing
    try 
        bond=_parse_bond_section(universe.sects["Bonds"],1)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse bondtype section."
    end
    bond[1]
end

function angles(universe)
    bond=nothing
    try 
        bond=_parse_bond_section(universe.sects["Bonds"],2)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse angles section."
    end
    bond[2]
end

function angletype(universe)
    bond=nothing
    try 
        bond=_parse_bond_section(universe.sects["Bonds"],2)
    catch e
        @error "LAMMPSDataReader.jl Failed to parse angletype section."
    end
    bond[1]
end

#uni=Universe("./example/input2.data")
#print(box(uni))
#print(charges(uni))