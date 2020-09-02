<?php

#################
### FUNCTIONS ###
#################



### function pp
###############

function pp($arr){
    $retStr = '<ul>';
    if (is_array($arr)){
        foreach ($arr as $key=>$val){
            if (is_array($val)){
                $retStr .= '<li>' . $key . ' => ' . pp($val) . '</li>';
            }else{
                $retStr .= '<li>' . $key . ' => ' . $val . '</li>';
            }
        }
    }
    $retStr .= '</ul>';
    return $retStr;
}



### function path_full
######################

function path_full($PATH="") {
    $PATH_SPLIT=explode ( "/" , $PATH );
    if (isset($PATH_SPLIT[0]) && $PATH_SPLIT[0]!="") {$ROOT=$PATH_SPLIT[0];} else {$ROOT="*";};
    if (isset($PATH_SPLIT[1]) && $PATH_SPLIT[1]!="") {$REPOSITORY=$PATH_SPLIT[1];} else {$REPOSITORY="*";};
    if (isset($PATH_SPLIT[2]) && $PATH_SPLIT[2]!="") {$GROUP=$PATH_SPLIT[2];} else {$GROUP="*";};
    if (isset($PATH_SPLIT[3]) && $PATH_SPLIT[3]!="") {$PROJECT=$PATH_SPLIT[3];} else {$PROJECT="*";};
    if (isset($PATH_SPLIT[4]) && $PATH_SPLIT[4]!="") {$RUN=$PATH_SPLIT[4];} else {$RUN="*";};
    if (isset($PATH_SPLIT[5]) && $PATH_SPLIT[5]!="") {$SAMPLE=$PATH_SPLIT[5];} else {$SAMPLE="*";};
    return "$ROOT/$REPOSITORY/$GROUP/$PROJECT/$RUN/$SAMPLE/";
};



### function path_short
#######################

function path_short($PATH="") {
    return str_replace("//","/",str_replace("*/","",path_full($PATH)));
};



### function path_html
######################

function path_html($PATH="") {
    # <a href='?PATH=$REPOSITORY/$GROUP/$PROJECT/$RUN/$SAMPLE/'>$SAMPLE</a>
    $PATH_SHORT=path_short($PATH);
    $PATH_HREF="";

    #$last_value="";
    foreach (explode("/",$PATH_SHORT) as $key => $value) {

        if ($value!="") {

            $PATH_HREF.="$value/";
            if ($key==0) {
                $PATH_RETURN="
                    <big>
                    <a class='' href='?PATH=$PATH_HREF'>
                        <div class='card-img align-self-center mbr-bold'>
                            <span class='mbr-iconfont mbri-home mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
                        </div>
                        HOME
                    </a>
                    </big>
                ";
            } else {

                #$PATH_RETURN.=" > <a class='btn btn-primary-outline display-7' href='?PATH=$PATH_HREF'>$value</a>";
                $PATH_RETURN.="
                &nbsp;
                &nbsp;
                <a class='' href='?PATH=$PATH_HREF'>
                    <div class='align-self-center bold'>
                        <span class='mbr-iconfont mbri-arrow-next mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
                        <span class='mbr-iconfont mbri-folder mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>
                    </div>
                    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$value
                </a>
                &nbsp;
                ";
            };
            # card-img
        };
        #$last_value=$value;
    };

    if ($key<6) {
        $PATH_RETURN.="

        &nbsp;
        &nbsp;
        <a class='' href='?PATH=$PATH_HREF'>
            <div class='align-self-center bold'>

                <span class='mbr-iconfont mbri-arrow-down mbr-bold' style='color: rgb(20, 157, 204); fill: rgb(20, 157, 204);'></span>

                &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;".path_links($PATH)."
            </div>

        </a>
        &nbsp;

        ";
    }

    # mbri-cust-feedback

    return $PATH_RETURN;
};



### function path_links
#######################

function path_links($PATH="") {
    $PATH_SHORT=path_short($PATH);
    if (count(explode( "/", $PATH ))<=6) {
        foreach (glob ( $PATH_SHORT."/*", GLOB_ONLYDIR ) as $key => $value) {
            $PATH_LINKS.="
                <a class='' href='?PATH="."$PATH_SHORT".end(explode( "/", $value ))."/'>
                    <div class=' align-self-left align-left bold'>
                        ".end(explode( "/", $value ))."
                    </div>
                </a>
            ";
        };
    };
    #echo "L $LINKS_HTML L"; btn btn-primary-outline display-8 display-8
    return $PATH_LINKS;
};



### function rutime
###################

function rutime($ru, $rus, $index) {
    return ($ru["ru_$index.tv_sec"]*1000 + intval($ru["ru_$index.tv_usec"]/1000))
     -  ($rus["ru_$index.tv_sec"]*1000 + intval($rus["ru_$index.tv_usec"]/1000));
}



### function tags_extract
#########################

function tags_extract($TAGS="") {
    $return=str_replace(" ","\n",str_replace("!"," ",$TAGS));
    preg_match_all('/.*#.*/i', $return, $matchWords);
    return implode(" ",$matchWords[0]);
    #  $matchWords = $matchWords[0];
    #  return $return;
};







?>