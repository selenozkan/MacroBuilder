header_html = '''
<!DOCTYPE html>
<html lang="en">

<head>
  <meta charset="utf-8">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
  <title>Superimposition Report</title>
  <!-- Bootstrap core CSS-->
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/css/bootstrap.min.css">
  <link href="css/sb-admin.css" rel="stylesheet">
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.3.7/css/bootstrap.min.css" />
  <link rel="stylesheet" href="https://drvic10k.github.io/bootstrap-sortable/Contents/bootstrap-sortable.css" />
  <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.2.1/jquery.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.3.7/js/bootstrap.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/moment.js/2.19.1/moment.js"></script>
  <script src="https://drvic10k.github.io/bootstrap-sortable/Scripts/bootstrap-sortable.js"></script>
  <link href="//netdna.bootstrapcdn.com/font-awesome/4.0.3/css/font-awesome.css" rel="stylesheet">

<script> {ngl} </script>

  <style>
  div.solid {{border-style: solid; border-color:#bababa;}}
  html {{position: relative;min-height: 100%;}}
  body {{overflow-x: hidden; font-size: 10pt;line-height: 25px; }}
  div.content-wrapper {{min-height: calc(100vh - 56px);padding-top: 1rem;overflow-x: hidden;background: white;}}
  #myBtn {{display: none;position: fixed;bottom: 20px;right: 20px;z-index: 99;font-size: 15px;border: none;outline: none;
    background-color: #5342f4;color: white;cursor: pointer;padding: 15px;border-radius: 3px;}}
  #myBtn:hover {{background-color: #555;}}
  div.card-header{{font-weight: bold;}}
  .mol-container {{width: 100%;height: 300px;position: relative;}}
  .card-body-plot {{text-align: center}}
</style>

</head>
<body id="page-top">
  <div class="row">
    <div class="col-sm-1"></div>
    <div class="col-sm-10">
  <button onclick="topFunction()" id="myBtn" title="Go to top">Top</button>
  <div class="content-wrapper">
    <div class="container-fluid">
      <!-- Breadcrumbs-->
      <ol class="breadcrumb">
        <li class="breadcrumb-item active"><h3>SUPERIMPOSITION by PROTEIN STRUCTURE</h3></li>
      </ol>
      <!-- Job-->
'''

job_html = '''
<div class="card mb-3">
<div class="card-header">Job</div>
<div class="card-body">
<ul>
<li><b>Commandline arguments:</b> {arguments} </li>
<li><b>Inputfile:</b>
 <a class="btn btn-xs btn-default" data-toggle="collapse" data-target="#demo"> <span class="glyphicon glyphicon-eye-open"></span>
	</a>
<div id="demo" class="collapse">
<pre>
{pairs_file}
</pre>
</div>

<ul>
<li><b>number of pairs introduced:</b> {n_pairs}</li>
</ul>
</li>
</ul>
</div>
</div>
<ol class="breadcrumb">
<li class="breadcrumb-item active"><h5>RESULT</h5></li>
</ol>
'''
table_header = """
<table class="table table-striped sortable" >
        <thead>
            <tr>
                <th>Intermediate name</th>
                <th>Pairs superimposed</th>
                <th>RMSD</th>
                <th>min distance*</th>
                <th>max distance*</th>
                <th></th>
            </tr>
        </thead>
        <tfoot>
            <tr>
            <th>Intermediate name</th>
            <th>Pairs superimposed</th>
            <th>RMSD</th>
            <th>min distance</th>
            <th>max distance</th>
            <th></th>
            </tr>
        </tfoot>
        <tbody>
    """


table_entry ="""
            <tr>
                <td> {complex_name} </td>
                <td> {pair_1} , {pair_2}</td>
                <td>{RMSD}</td>
                <td>{min_dist}</td>
                <td>{max_dist} </td>
                <td><a href = {complex_name}.html target=_black> <span class="glyphicon glyphicon-eye-open"></span></a><td>
            </tr>
        """

table_bot= """
    </tbody>
    </table>
    <small>* Distance between previous intermediate complex (base) and the added structure (not the superimposed residues, those are not part of the new complex).</small>
         """


final_html='''
      <div class="card mb-2">
        <div class="card-header">Final Complex</div>
        <div class="card-body">
          <div class="row">
            <div class="col-sm-2"></div>
            <div class="col-sm-8">
              <script>
                  document.addEventListener("DOMContentLoaded", function () {{
                    var stage = new NGL.Stage("viewport_final");
                    stage.loadFile("{structure_file}", {{defaultRepresentation: true}});
                  }});
              </script>
                <div id="viewport_final" style="width:650px; height:500px;"></div>
              </div>
              <div class="col-sm-2"></div>
      </div>
      </div>
    </div>
      <!-- Area Chart Example-->
      <div class="card mb-3">
        <div class="card-header">
          <i class="fa fa-area-chart"></i>Protein Distance Plot</div>
        <div class="card-body-plot">
          <img src='{plot_file}' width="90%" height="100%"s>
        </div>
      </div>
      <ol class="breadcrumb">
      <li class="breadcrumb-item active"><h5>INTERMEDIATE COMPLEXES:</h5></li>
      </ol>
'''

base_html ='''
      <br>
      <ol class="breadcrumb">
      <li class="breadcrumb-item active"><h5>INTERMEDIATE COMPLEX:</h5></li>
      </ol>
      <div class="bases">
      <div class="solid">
      <div class="card mb-3">
        <div class="card-header">Complex name: {complex_name}</div>
        <div class="card-body">
          <div class="row">
            <div class="col-sm-6" style="padding-left:40px">
                Superimposition completed successfully.
                <ul>
                <li>Pairs being superimposed: {pairs}.</li>
                <li>Root mean square deviation (RMSD): {RMSD}</li>
                </ul>
                Between previous base and added residues:
                <ul>
                <li>Minimum distance (Amstrong): {min_dist} </li>
                <li>Maximum distance (Amstrong): {max_dist} </li>
              </ul>
            </div>
            <div class="col-sm-5">
              <script>
                  document.addEventListener("DOMContentLoaded", function () {{
                    var stage = new NGL.Stage("viewport_{n}");
                    var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {{
                    this.atomColor = function (atom) {{
                    if ({chain_names_base}.includes(atom.chainname)) {{
                    return 0xFFE5BC; // red
                    }}else {{
                    return 0xFF0202;  //
                    }}
                    }};
                    }});
        stage.loadFile("{structure_f}").then(function (o) {{
        o.addRepresentation("cartoon", {{ color: schemeId }})
        o.autoView();}})
                  }});

              </script>
                <div id="viewport_{n}" style="width:450px; height:400px;"></div>
            </div>
          </div>
      </div>
      </div>
      <!-- Area Chart Example-->
      <div class="card mb-3">
        <div class="card-header">
          <i class="fa fa-area-chart"></i> Protein Distance Plot</div>
        <div class="card-body-plot">
          <img src='{plot_f}' width="100%">
        </div>
      </div>
      </div>'''

footer_html='''
      <br>
      <br>

    <!-- /.container-fluid-->
    <footer class="sticky-footer">
      <div class="container">
        <div class="text-center">
          <small>Version 1.0 / Last update: March,2018 / Antía Fernández Pintos, Eva Martín del Pico, Selen Özkan<br></small>
        </div>
        <div class="text-center">
          <small>Universitat Pompeu Fabra, Bioinformatics For Health Sciences, 2018<br></small>
        </div>
      </div>
    </footer>

    <script>
    // When the user scrolls down 20px from the top of the document, show the button
    window.onscroll = function() {scrollFunction()};

    function scrollFunction() {
        if (document.body.scrollTop > 20 || document.documentElement.scrollTop > 20) {
            document.getElementById("myBtn").style.display = "block";
        } else {
            document.getElementById("myBtn").style.display = "none";
        }
    }

    // When the user clicks on the button, scroll to the top of the document
    function topFunction() {
        document.body.scrollTop = 0;
        document.documentElement.scrollTop = 0;
    }
    </script>

  </div>
  <!-- /.content-wrapper-->
</div>
<div class="col-sm-2"></div>
</div>
</body>

</html>'''
