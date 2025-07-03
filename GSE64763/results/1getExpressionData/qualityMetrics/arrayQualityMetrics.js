// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "GSM1579399_DL_U133A_2_FBR001.CEL.gz", "GSM1579399_DL_U133A_2_FBR001.CEL.gz", "Fibroid", "Fibroid", "10/26/06 12:17:08" ], [ "2", "GSM1579400_DL_U133A_2_FBR002.CEL.gz", "GSM1579400_DL_U133A_2_FBR002.CEL.gz", "Fibroid", "Fibroid", "10/26/06 12:32:24" ], [ "3", "GSM1579401_DL_U133A_2_FBR005.CEL.gz", "GSM1579401_DL_U133A_2_FBR005.CEL.gz", "Fibroid", "Fibroid", "10/26/06 11:54:23" ], [ "4", "GSM1579402_DL_U133A_2_FBR006.CEL.gz", "GSM1579402_DL_U133A_2_FBR006.CEL.gz", "Fibroid", "Fibroid", "10/26/06 11:46:47" ], [ "5", "GSM1579403_DL_U133A_2_FBR007.CEL.gz", "GSM1579403_DL_U133A_2_FBR007.CEL.gz", "Fibroid", "Fibroid", "10/26/06 12:01:59" ], [ "6", "GSM1579404_DL_U133A_2_FBR008.CEL.gz", "GSM1579404_DL_U133A_2_FBR008.CEL.gz", "Fibroid", "Fibroid", "10/26/06 12:50:33" ], [ "7", "GSM1579405_DL_U133A_2_FBR009.CEL.gz", "GSM1579405_DL_U133A_2_FBR009.CEL.gz", "Fibroid", "Fibroid", "10/26/06 12:58:12" ], [ "8", "GSM1579406_DL_U133A_2_FBR0010.CEL.gz", "GSM1579406_DL_U133A_2_FBR0010.CEL.gz", "Fibroid", "Fibroid", "10/26/06 12:24:47" ], [ "9", "GSM1579407_DL_U133A_2_FBR011.CEL.gz", "GSM1579407_DL_U133A_2_FBR011.CEL.gz", "Fibroid", "Fibroid", "10/26/06 13:05:57" ], [ "10", "GSM1579408_DL_U133A_2_FBR013.CEL.gz", "GSM1579408_DL_U133A_2_FBR013.CEL.gz", "Fibroid", "Fibroid", "10/26/06 13:14:17" ], [ "11", "GSM1579409_DL_U133A_2_FBR015.CEL.gz", "GSM1579409_DL_U133A_2_FBR015.CEL.gz", "Fibroid", "Fibroid", "10/26/06 12:42:59" ], [ "12", "GSM1579410_DL_U133A_2_FBR016.CEL.gz", "GSM1579410_DL_U133A_2_FBR016.CEL.gz", "Fibroid", "Fibroid", "10/25/06 13:03:32" ], [ "13", "GSM1579411_DL_U133A_2_FBR019.CEL.gz", "GSM1579411_DL_U133A_2_FBR019.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:32:48" ], [ "14", "GSM1579412_DL_U133A_2_FBR020.CEL.gz", "GSM1579412_DL_U133A_2_FBR020.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:55:49" ], [ "15", "GSM1579413_DL_U133A_2_FBR021.CEL.gz", "GSM1579413_DL_U133A_2_FBR021.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:17:12" ], [ "16", "GSM1579414_DL_U133A_2_FBR022.CEL.gz", "GSM1579414_DL_U133A_2_FBR022.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:24:54" ], [ "17", "GSM1579415_DL_U133A_2_FBR024.CEL.gz", "GSM1579415_DL_U133A_2_FBR024.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:09:27" ], [ "18", "GSM1579416_DL_U133A_2_FBR025.CEL.gz", "GSM1579416_DL_U133A_2_FBR025.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:01:39" ], [ "19", "GSM1579417_DL_U133A_2_FBR027.CEL.gz", "GSM1579417_DL_U133A_2_FBR027.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:48:04" ], [ "20", "GSM1579418_DL_U133A_2_FBR028.CEL.gz", "GSM1579418_DL_U133A_2_FBR028.CEL.gz", "Fibroid", "Fibroid", "10/25/06 12:40:20" ], [ "21", "GSM1579419_DL_U133A_2_fbr_030.CEL.gz", "GSM1579419_DL_U133A_2_fbr_030.CEL.gz", "Fibroid", "Fibroid", "02/28/07 13:02:56" ], [ "22", "GSM1579420_DL_U133A_2_fbr_031.CEL.gz", "GSM1579420_DL_U133A_2_fbr_031.CEL.gz", "Fibroid", "Fibroid", "02/28/07 13:10:38" ], [ "23", "GSM1579421_DL_U133A_2_fbr_033.CEL.gz", "GSM1579421_DL_U133A_2_fbr_033.CEL.gz", "Fibroid", "Fibroid", "02/28/07 13:18:26" ], [ "24", "GSM1579422_DL_U133A_2_fbr_034.CEL.gz", "GSM1579422_DL_U133A_2_fbr_034.CEL.gz", "Fibroid", "Fibroid", "02/28/07 13:26:12" ], [ "25", "GSM1579423_DL_U133A_2_fbr_036.CEL.gz", "GSM1579423_DL_U133A_2_fbr_036.CEL.gz", "Fibroid", "Fibroid", "02/28/07 13:33:54" ], [ "26", "GSM1579424_DL_U133A_2_JB1075.CEL.gz", "GSM1579424_DL_U133A_2_JB1075.CEL.gz", "Tumor", "Tumor", "10/25/06 11:53:49" ], [ "27", "GSM1579425_DL_U133A_2_JB109.CEL.gz", "GSM1579425_DL_U133A_2_JB109.CEL.gz", "Tumor", "Tumor", "10/26/06 11:39:01" ], [ "28", "GSM1579426_DL_U133A_2_JB1541.CEL.gz", "GSM1579426_DL_U133A_2_JB1541.CEL.gz", "Tumor", "Tumor", "10/26/06 12:09:34" ], [ "29", "GSM1579427_DL_U133A_2_JB1583.CEL.gz", "GSM1579427_DL_U133A_2_JB1583.CEL.gz", "Tumor", "Tumor", "11/03/06 12:39:31" ], [ "30", "GSM1579428_DL_U133A_2_jb_2124.CEL.gz", "GSM1579428_DL_U133A_2_jb_2124.CEL.gz", "Tumor", "Tumor", "02/28/07 13:41:34" ], [ "31", "GSM1579429_DL_U133A_2_JB2486.CEL.gz", "GSM1579429_DL_U133A_2_JB2486.CEL.gz", "Tumor", "Tumor", "11/03/06 12:58:47" ], [ "32", "GSM1579430_DL_U133A_2_JB2692.CEL.gz", "GSM1579430_DL_U133A_2_JB2692.CEL.gz", "Tumor", "Tumor", "11/03/06 13:06:26" ], [ "33", "GSM1579431_DL_U133A_2_jb_2883.CEL.gz", "GSM1579431_DL_U133A_2_jb_2883.CEL.gz", "Tumor", "Tumor", "02/28/07 13:49:29" ], [ "34", "GSM1579432_DL_U133A_2_JB72.CEL.gz", "GSM1579432_DL_U133A_2_JB72.CEL.gz", "Tumor", "Tumor", "11/03/06 13:22:08" ], [ "35", "GSM1579433_DL_U133A_2_JB804.CEL.gz", "GSM1579433_DL_U133A_2_JB804.CEL.gz", "Tumor", "Tumor", "11/03/06 13:14:31" ], [ "36", "GSM1579434_DL_U133A_2_tb_117247.CEL.gz", "GSM1579434_DL_U133A_2_tb_117247.CEL.gz", "Tumor", "Tumor", "02/28/07 15:03:21" ], [ "37", "GSM1579435_DL_U133A_2_tb_120946.CEL.gz", "GSM1579435_DL_U133A_2_tb_120946.CEL.gz", "Tumor", "Tumor", "02/28/07 14:55:32" ], [ "38", "GSM1579436_DL_U133A_2_tb_14004.CEL.gz", "GSM1579436_DL_U133A_2_tb_14004.CEL.gz", "Tumor", "Tumor", "02/28/07 15:18:49" ], [ "39", "GSM1579437_DL_U133A_2_tb_6424.CEL.gz", "GSM1579437_DL_U133A_2_tb_6424.CEL.gz", "Tumor", "Tumor", "02/28/07 15:26:31" ], [ "40", "GSM1579438_DL_U133A_2_tb_9004.CEL.gz", "GSM1579438_DL_U133A_2_tb_9004.CEL.gz", "Tumor", "Tumor", "02/28/07 15:34:09" ], [ "41", "GSM1579439_DL_U133A_2_32.CEL.gz", "GSM1579439_DL_U133A_2_32.CEL.gz", "Tumor", "Tumor", "02/02/07 11:37:50" ], [ "42", "GSM1579440_DL_U133A_2_34.CEL.gz", "GSM1579440_DL_U133A_2_34.CEL.gz", "Tumor", "Tumor", "02/02/07 11:29:58" ], [ "43", "GSM1579441_DL_U133A_2_36.CEL.gz", "GSM1579441_DL_U133A_2_36.CEL.gz", "Tumor", "Tumor", "02/02/07 11:45:40" ], [ "44", "GSM1579442_DL_U133A_2_40.CEL.gz", "GSM1579442_DL_U133A_2_40.CEL.gz", "Tumor", "Tumor", "02/02/07 11:22:17" ], [ "45", "GSM1579443_DL_U133A_2_71.CEL.gz", "GSM1579443_DL_U133A_2_71.CEL.gz", "Tumor", "Tumor", "11/15/06 12:39:52" ], [ "46", "GSM1579444_DL_U133A_2_72.CEL.gz", "GSM1579444_DL_U133A_2_72.CEL.gz", "Tumor", "Tumor", "11/15/06 11:53:55" ], [ "47", "GSM1579445_DL_U133A_2_73.CEL.gz", "GSM1579445_DL_U133A_2_73.CEL.gz", "Tumor", "Tumor", "11/15/06 12:32:14" ], [ "48", "GSM1579446_DL_U133A_2_74.CEL.gz", "GSM1579446_DL_U133A_2_74.CEL.gz", "Tumor", "Tumor", "11/15/06 12:01:35" ], [ "49", "GSM1579447_DL_U133A_2_75.CEL.gz", "GSM1579447_DL_U133A_2_75.CEL.gz", "Tumor", "Tumor", "11/15/06 13:02:51" ], [ "50", "GSM1579448_DL_U133A_2_76.CEL.gz", "GSM1579448_DL_U133A_2_76.CEL.gz", "Tumor", "Tumor", "11/15/06 13:18:03" ], [ "51", "GSM1579449_DL_U133A_2_NL001.CEL.gz", "GSM1579449_DL_U133A_2_NL001.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 11:51:33" ], [ "52", "GSM1579450_DL_U133A_2_NL002.CEL.gz", "GSM1579450_DL_U133A_2_NL002.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 12:15:10" ], [ "53", "GSM1579451_DL_U133A_2_NL003.CEL.gz", "GSM1579451_DL_U133A_2_NL003.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 11:59:26" ], [ "54", "GSM1579452_DL_U133A_2_NL004.CEL.gz", "GSM1579452_DL_U133A_2_NL004.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 12:23:13" ], [ "55", "GSM1579453_DL_U133A_2_NL005.CEL.gz", "GSM1579453_DL_U133A_2_NL005.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 12:31:20" ], [ "56", "GSM1579454_DL_U133A_2_NL006.CEL.gz", "GSM1579454_DL_U133A_2_NL006.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 11:35:32" ], [ "57", "GSM1579455_DL_U133A_2_NL007.CEL.gz", "GSM1579455_DL_U133A_2_NL007.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 11:43:40" ], [ "58", "GSM1579456_DL_U133A_2_NL009.CEL.gz", "GSM1579456_DL_U133A_2_NL009.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 11:26:07" ], [ "59", "GSM1579457_DL_U133A_2_NL010.CEL.gz", "GSM1579457_DL_U133A_2_NL010.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 11:18:12" ], [ "60", "GSM1579458_DL_U133A_2_NL011.CEL.gz", "GSM1579458_DL_U133A_2_NL011.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 12:07:24" ], [ "61", "GSM1579459_DL_U133A_2_NL012.CEL.gz", "GSM1579459_DL_U133A_2_NL012.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/03/06 12:47:37" ], [ "62", "GSM1579460_DL_U133A_2_NL013.CEL.gz", "GSM1579460_DL_U133A_2_NL013.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 11:51:00" ], [ "63", "GSM1579461_DL_U133A_2_NL015.CEL.gz", "GSM1579461_DL_U133A_2_NL015.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 11:58:34" ], [ "64", "GSM1579462_DL_U133A_2_NL018.CEL.gz", "GSM1579462_DL_U133A_2_NL018.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 12:06:20" ], [ "65", "GSM1579463_DL_U133A_2_NL019.CEL.gz", "GSM1579463_DL_U133A_2_NL019.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 12:29:47" ], [ "66", "GSM1579464_DL_U133A_2_NL021.CEL.gz", "GSM1579464_DL_U133A_2_NL021.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 12:21:55" ], [ "67", "GSM1579465_DL_U133A_2_NL022.CEL.gz", "GSM1579465_DL_U133A_2_NL022.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 12:14:07" ], [ "68", "GSM1579466_DL_U133A_2_NL023.CEL.gz", "GSM1579466_DL_U133A_2_NL023.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 12:37:24" ], [ "69", "GSM1579467_DL_U133A_2_NL024.CEL.gz", "GSM1579467_DL_U133A_2_NL024.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/10/06 12:45:05" ], [ "70", "GSM1579468_DL_U133A_2_NL025.CEL.gz", "GSM1579468_DL_U133A_2_NL025.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/15/06 12:24:34" ], [ "71", "GSM1579469_DL_U133A_2_NL026.CEL.gz", "GSM1579469_DL_U133A_2_NL026.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/15/06 12:16:52" ], [ "72", "GSM1579470_DL_U133A_2_NL027.CEL.gz", "GSM1579470_DL_U133A_2_NL027.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/15/06 12:09:17" ], [ "73", "GSM1579471_DL_U133A_2_NL028.CEL.gz", "GSM1579471_DL_U133A_2_NL028.CEL.gz", "NormalMyometrium", "NormalMyometrium", "11/15/06 12:47:31" ], [ "74", "GSM1579472_DL_U133A_2_nl_029.CEL.gz", "GSM1579472_DL_U133A_2_nl_029.CEL.gz", "NormalMyometrium", "NormalMyometrium", "02/28/07 13:57:09" ], [ "75", "GSM1579473_DL_U133A_2_nl_030.CEL.gz", "GSM1579473_DL_U133A_2_nl_030.CEL.gz", "NormalMyometrium", "NormalMyometrium", "02/28/07 14:04:54" ], [ "76", "GSM1579474_DL_U133A_2_nl_031.CEL.gz", "GSM1579474_DL_U133A_2_nl_031.CEL.gz", "NormalMyometrium", "NormalMyometrium", "02/28/07 14:12:38" ], [ "77", "GSM1579475_DL_U133A_2_nl_033.CEL.gz", "GSM1579475_DL_U133A_2_nl_033.CEL.gz", "NormalMyometrium", "NormalMyometrium", "02/28/07 14:20:25" ], [ "78", "GSM1579476_DL_U133A_2_nl_035.CEL.gz", "GSM1579476_DL_U133A_2_nl_035.CEL.gz", "NormalMyometrium", "NormalMyometrium", "02/28/07 14:28:07" ], [ "79", "GSM1579477_DL_U133A_2_nl_036.CEL.gz", "GSM1579477_DL_U133A_2_nl_036.CEL.gz", "NormalMyometrium", "NormalMyometrium", "02/28/07 15:11:05" ] ];
var svgObjectNames   = [ "pca", "dens", "dig" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
