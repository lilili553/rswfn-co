pub static GJF: &str = r#"%chk={chk}.chk
# {job}

Generate by Pywfn

{ele.0} {ele.1}
{geome}


"#;

pub static XYZ: &str = r#"{natm}
{title}
{geome}
"#;

pub static SI: &str = r#"{name}
{engs}

Cartesian coordinate
{geome}

Variable frequencies
{freqs}
"#;
