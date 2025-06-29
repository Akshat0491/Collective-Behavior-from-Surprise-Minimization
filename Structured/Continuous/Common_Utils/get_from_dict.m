function ret=get_from_dict(dict1,var_name) %dict maps string to a cell containing deserived object. "Name"-->{Name}
    dummy=dict1(var_name);
    ret=dummy{1};
end