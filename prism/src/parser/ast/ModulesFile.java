//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package parser.ast;

import java.util.*;

import parser.*;
import parser.visitor.*;
import prism.PrismLangException;
import prism.ModelType;
import parser.type.*;

// Class representing parsed model file

public class ModulesFile extends ASTElement
{
	// Model type (enum)
	private ModelType modelType;

	// Model components
	private FormulaList formulaList;
	private LabelList labelList;
	private ConstantList constantList;
	private Vector<Declaration> globals; // Global variables
	private Vector<Object> modules; // Modules (includes renamed modules)
	private SystemDefn systemDefn; // System definition (system...endsystem construct)
	private ArrayList<RewardStruct> rewardStructs; // Rewards structures
	private Expression initStates; // Initial states specification

	// Lists of all identifiers used
	private Vector<String> formulaIdents;
	private Vector<String> constantIdents;
	private Vector<String> varIdents; // TODO: don't need?
	// List of all module names
	private String[] moduleNames;
	// List of synchronising actions
	private Vector<String> synchs;
	//TODO quick workaround for FDCTMC models- when synchronization will be done properly for 
	// FDCTMC models, the variable below probably won't be needed anymore
	//mapping of all actions to modules where they occur
	private Map<String,Integer> synchsOfFDToModule;
	private Map<String,Integer> synchsOfFDToFD;
	// Lists of variable info (declaration, name, type)
	private Vector<Declaration> varDecls;
	private Vector<String> varNames;
	private Vector<Type> varTypes;

	// actual values of (some or all) constants
	private Values constantValues;

	// Constructor

	public ModulesFile()
	{
		formulaList = new FormulaList();
		labelList = new LabelList();
		constantList = new ConstantList();
		modelType = ModelType.MDP; // default type
		globals = new Vector<Declaration>();
		modules = new Vector<Object>();
		systemDefn = null;
		rewardStructs = new ArrayList<RewardStruct>();
		initStates = null;
		formulaIdents = new Vector<String>();
		constantIdents = new Vector<String>();
		varIdents = new Vector<String>();
		varDecls = new Vector<Declaration>();
		varNames = new Vector<String>();
		varTypes = new Vector<Type>();
		constantValues = null;
	}

	// Set methods

	public void setFormulaList(FormulaList fl)
	{
		formulaList = fl;
	}

	public void setLabelList(LabelList ll)
	{
		labelList = ll;
	}

	public void setConstantList(ConstantList cl)
	{
		constantList = cl;
	}

	public void setModelType(ModelType t)
	{
		modelType = t;
	}

	public void addGlobal(Declaration d)
	{
		globals.add(d);
	}

	public void setGlobal(int i, Declaration d)
	{
		globals.set(i, d);
	}

	public void addModule(Module m)
	{
		modules.add(m);
		m.setParent(this);
	}

	public void setModule(int i, Module m)
	{
		modules.set(i, m);
		m.setParent(this);
	}

	public void addRenamedModule(RenamedModule m)
	{
		modules.add(m);
	}

	public void setSystemDefn(SystemDefn s)
	{
		systemDefn = s;
	}

	public void clearRewardStructs()
	{
		rewardStructs.clear();
	}

	public void addRewardStruct(RewardStruct r)
	{
		rewardStructs.add(r);
	}

	public void setRewardStruct(int i, RewardStruct r)
	{
		rewardStructs.set(i, r);
	}

	// this method is included for backwards compatibility only
	public void setRewardStruct(RewardStruct r)
	{
		rewardStructs.clear();
		rewardStructs.add(r);
	}

	public void setInitialStates(Expression e)
	{
		initStates = e;
	}

	// Get methods

	public FormulaList getFormulaList()
	{
		return formulaList;
	}

	public LabelList getLabelList()
	{
		return labelList;
	}

	public ConstantList getConstantList()
	{
		return constantList;
	}

	public ModelType getModelType()
	{
		return modelType;
	}

	public String getTypeString()
	{
		return "" + modelType;
	}

	public String getTypeFullString()
	{
		return modelType.fullName();
	}

	public int getNumGlobals()
	{
		return globals.size();
	}

	public Declaration getGlobal(int i)
	{
		return globals.elementAt(i);
	}

	public int getNumModules()
	{
		return modules.size();
	}

	// get module by index
	// returns null if it is a RenamedModule
	// these will have been replaced with Modules after tidyUp()
	public Module getModule(int i)
	{
		Object o = modules.elementAt(i);
		return (o instanceof Module) ? (Module) o : null;
	}

	// get the index of a module by its name
	// (returns -1 if it does not exist)
	// (or is a RenamedModule and hasn't been turned into a normal Module yet)
	public int getModuleIndex(String s)
	{
		int i;
		Module m;

		for (i = 0; i < modules.size(); i++) {
			m = getModule(i);
			if (m != null)
				if (s.equals(m.getName())) {
					return i;
				}
		}

		return -1;
	}

	public SystemDefn getSystemDefn()
	{
		return systemDefn;
	}

	/**
	 * Get the number of reward structures in the model.
	 */
	public int getNumRewardStructs()
	{
		return rewardStructs.size();
	}

	/**
	 * Get a reward structure by its index
	 * (indexed from 0, not from 1 like at the user (property language) level).
	 * Returns null if index is out of range.
	 */
	public RewardStruct getRewardStruct(int i)
	{
		return (i < rewardStructs.size()) ? rewardStructs.get(i) : null;
	}

	/**
	 * Get access to the list of reward structures
	 */
	public List<RewardStruct> getRewardStructs()
	{
		return rewardStructs;
	}

	/**
	 * Get the index of a module by its name
	 * (indexed from 0, not from 1 like at the user (property language) level).
	 * Returns -1 if name does not exist.
	 */
	public int getRewardStructIndex(String name)
	{
		int i, n;
		n = rewardStructs.size();
		for (i = 0; i < n; i++) {
			if ((rewardStructs.get(i)).getName().equals(name))
				return i;
		}
		return -1;
	}

	/**
	 * Get a reward structure by its name
	 * Returns null if name does not exist.
	 */
	public RewardStruct getRewardStructByName(String name)
	{
		int i = getRewardStructIndex(name);
		return i == -1 ? null : getRewardStruct(i);
	}

	/**
	 * Get the first reward structure (exists for backwards compatibility only).
	 */
	public RewardStruct getRewardStruct()
	{
		return getRewardStruct(0);
	}

	/**
	 * Get the expression used in init...endinit to define the initial states of the model.
	 * This is null if absent, i.e. the model has a single initial state defined by the
	 * initial values for each individual variable.  
	 * If non-null, we have to assume that there may be multiple initial states. 
	 */
	public Expression getInitialStates()
	{
		return initStates;
	}

	/**
	 * Look up a property by name.
	 * Returns null if not found.
	 * Currently only exists for forwards compatibility.
	 */
	public Property getPropertyByName(String name)
	{
		return null;
	}
	
	/**
	 * Check if an identifier is used by this model
	 * (as a formula, constant, or variable)
	 */
	public boolean isIdentUsed(String ident)
	{
		return formulaIdents.contains(ident) || constantIdents.contains(ident) || varIdents.contains(ident);
	}

	// get individual module name
	public String getModuleName(int i)
	{
		return moduleNames[i];
	}

	// get array of all module names
	public String[] getModuleNames()
	{
		return moduleNames;
	}

	/**
	 * Get the list of action names.
	 */
	public Vector<String> getSynchs()
	{
		return synchs;
	}

	/**
	 * Get the {@code i}th action name (0-indexed).
	 */
	public String getSynch(int i)
	{
		return synchs.get(i);
	}

	public boolean isSynch(String s)
	{
		if (synchs == null)
			return false;
		else
			return synchs.contains(s);
	}

	/**
	 * Get the any module index that contains given action label.
	 * @param label of synchronization action.
	 */	
	public int getModuleIndexFromFDSynch(String label){
		if(synchsOfFDToFD.get(label) == null) return 0;
		
		return synchsOfFDToModule.get(label);
	}
	
	/**
	 * Get the FD index that contains given action label.
	 * @param label of synchronization action.
	 */	
	public int getFDIndexFromSynch(String label){
		if(synchsOfFDToFD.get(label) == null) return -1;
		
		return synchsOfFDToFD.get(label);
	}
	
	// Variable query methods

	/**
	 * Get the total number of variables (global and local).
	 */
	public int getNumVars()
	{
		return varNames.size();
	}

	/**
	 * Look up the index of a variable in the model by name.
	 * Returns -1 if there is no such variable. 
	 */
	public int getVarIndex(String name)
	{
		return varNames.indexOf(name);
	}

	/**
	 * Get the declaration of the ith variable.
	 */
	public Declaration getVarDeclaration(int i)
	{
		return varDecls.get(i);
	}

	/**
	 * Get the name of the ith variable.
	 */
	public String getVarName(int i)
	{
		return varNames.get(i);
	}

	/**
	 * Get the type of the ith variable.
	 */
	public Type getVarType(int i)
	{
		return varTypes.get(i);
	}

	public Vector<String> getVarNames()
	{
		return varNames;
	}

	public Vector<Type> getVarTypes()
	{
		return varTypes;
	}

	public boolean isGlobalVariable(String s)
	{
		int i, n;

		n = getNumGlobals();
		for (i = 0; i < n; i++) {
			if (getGlobal(i).getName().equals(s))
				return true;
		}
		return false;
	}

	/**
	 * Method to "tidy up" after parsing (must be called)
	 * (do some checks and extract some information)
	 */
	public void tidyUp() throws PrismLangException
	{
		// Clear lists that will generated by this method 
		// (in case it has already been called previously).
		formulaIdents.clear();
		constantIdents.clear();
		varIdents.clear();
		varDecls.clear();
		varNames.clear();
		varTypes.clear();
		
		// Expansion of formulas and renaming

		// Check formula identifiers
		checkFormulaIdents();
		// Find all formulas (i.e. locate identifiers which are formulas).
		// Note: This should all be done before replacing any other identifiers
		// (e.g. with variables) because this relies on module renaming which in turn
		// must be done after formula expansion. Then, check for any cyclic
		// dependencies in the formula list and then expand all formulas.
		findAllFormulas(formulaList);
		formulaList.findCycles();
		expandFormulas(formulaList);
		// Perform module renaming
		sortRenamings();

		// Check label identifiers
		checkLabelIdents();

		// Check module names
		checkModuleNames();

		// Get synchronising action names
		getSynchNames();
		// Then identify/check any references to action names
		findAllActions(synchs);

		// Check constant identifiers
		checkConstantIdents();
		// Find all instances of constants
		// (i.e. locate identifiers which are constants)
		findAllConstants(constantList);
		// Check constants for cyclic dependencies
		constantList.findCycles();

		// Check variable names, etc.
		checkVarNames();
		// Find all instances of variables, replace identifiers with variables.
		// Also check variables valid, store indices, etc.
		findAllVars(varNames, varTypes);

		// Find all instances of property refs
		findAllPropRefs(this, null);
		
		// Check reward structure names
		checkRewardStructNames();

		// Various semantic checks 
		semanticCheck(this);
		// Type checking
		typeCheck();
		
		// If there are no undefined constants, set up values for constants
		// (to avoid need for a later call to setUndefinedConstants).
		// NB: Can't call setUndefinedConstants if there are undefined constants
		// because semanticCheckAfterConstants may fail. 
		if (getUndefinedConstants().isEmpty()) {
			setUndefinedConstants(null);
		}
	}

	// Check formula identifiers

	private void checkFormulaIdents() throws PrismLangException
	{
		int i, n;
		String s;

		n = formulaList.size();
		for (i = 0; i < n; i++) {
			s = formulaList.getFormulaName(i);
			if (isIdentUsed(s)) {
				throw new PrismLangException("Duplicated identifier \"" + s + "\"", formulaList.getFormulaNameIdent(i));
			} else {
				formulaIdents.add(s);
			}
		}
	}

	// sort out modules defined by renaming

	private void sortRenamings() throws PrismLangException
	{
		int i, j, n, i2, n2;
		RenamedModule module;
		Module newModule;
		String s;
		Object o;
		HashSet<String> renamedSoFar;

		// Go through modules and find ones which are defined by renaming
		n = modules.size();
		for (i = 0; i < n; i++) {
			o = modules.elementAt(i);
			if (o instanceof Module)
				continue;
			module = (RenamedModule) o;
			// Check base module exists
			j = getModuleIndex(module.getBaseModule());
			if (j == -1) {
				s = "No such module " + module.getBaseModule();
				s += " in renamed module \"" + module.getName() + "\"";
				throw new PrismLangException(s, module.getBaseModuleASTElement());
			}
			// Check for invalid renames
			n2 = module.getNumRenames();
			renamedSoFar = new HashSet<String>();
			for (i2 = 0; i2 < n2; i2++) {
				s = module.getOldName(i2);
				if (!renamedSoFar.add(s)) {
					throw new PrismLangException("Identifier \"" + s + "\" is renamed more than once in module \""
							+ module.getName() + "\"", module.getOldNameASTElement(i2));
				}
				if (formulaList.getFormulaIndex(s) != -1) {
					throw new PrismLangException("Formula \"" + s
							+ "\" cannot be renamed since formulas are expanded before module renaming", module
							.getOldNameASTElement(i2));
				}
			}
			// Then rename (a copy of) base module and replace
			// (note: also store name of base module for later reference)
			newModule = (Module) getModule(j).deepCopy().rename(module);
			newModule.setNameASTElement(module.getNameASTElement());
			newModule.setBaseModule(module.getBaseModule());
			setModule(i, newModule);
		}
	}

	// check label identifiers

	private void checkLabelIdents() throws PrismLangException
	{
		int i, n;
		String s;
		Vector<String> labelIdents;

		// go thru labels
		n = labelList.size();
		labelIdents = new Vector<String>();
		for (i = 0; i < n; i++) {
			s = labelList.getLabelName(i);
			// see if ident has been used already for a label
			if (labelIdents.contains(s)) {
				throw new PrismLangException("Duplicated label name \"" + s + "\"", labelList.getLabelNameIdent(i));
			} else {
				labelIdents.add(s);
			}
		}
	}

	// check module names

	private void checkModuleNames() throws PrismLangException
	{
		int i, j, n;
		String s;

		// check we have at least one module
		n = modules.size();
		if (n == 0) {
			throw new PrismLangException("There must be at least one module");
		}

		// compile list of all module names
		// and check as we go through
		moduleNames = new String[n];
		for (i = 0; i < n; i++) {
			s = getModule(i).getName();
			for (j = 0; j < i; j++) {
				if (s.equals(moduleNames[j])) {
					throw new PrismLangException("Duplicated module name \"" + s + "\"", getModule(i)
							.getNameASTElement());
				}
			}
			moduleNames[i] = s;
		}
	}

	// get all synch names

	private void getSynchNames() throws PrismLangException
	{
		Vector<String> v;
		String s;
		int i, j, n, m;

		// create vector to store names
		synchs = new Vector<String>();

		// go thru modules and extract names which appear in their commands
		n = modules.size();
		for (i = 0; i < n; i++) {
			v = getModule(i).getAllSynchs();
			m = v.size();
			for (j = 0; j < m; j++) {
				s = v.elementAt(j);
				if (!synchs.contains(s)) {
					synchs.add(s);
				}
			}
		}

		// then extract any which are introduced in system construct (i.e. by
		// renaming)
		if (systemDefn != null) {
			systemDefn.getSynchs(synchs);
		}
		
		if(!modelType.equals(prism.ModelType.FDCTMC)) return;
		
		synchsOfFDToModule = new HashMap<String, Integer>();
		synchsOfFDToFD = new HashMap<String, Integer>();
		for (i = 0; i < n; i++) {
			for (Command com : getModule(i).getCommands()) {
				String synchLabel = com.getSynch();
				int indexOfFD = ((Module)modules.get(i)).getFDelayIndex(com.getFDelay());
				if(synchLabel != null && indexOfFD != -1) {
					if(synchsOfFDToModule.get(synchLabel) != null && synchsOfFDToModule.get(synchLabel) != i)
						System.out.println("WARNING: Found multiple models synchronized on FD transition! Not supported yet!");
					if(synchsOfFDToFD.get(synchLabel) != null && synchsOfFDToFD.get(synchLabel) != indexOfFD)
						System.out.println("WARNING: Found multiple synchronized FD! Not supported yet!");
					synchsOfFDToModule.put(synchLabel, i);
					synchsOfFDToFD.put(synchLabel, indexOfFD);
				}
				
			}
		}
	}

	// check constant identifiers

	private void checkConstantIdents() throws PrismLangException
	{
		int i, n;
		String s;

		n = constantList.size();
		for (i = 0; i < n; i++) {
			s = constantList.getConstantName(i);
			if (isIdentUsed(s)) {
				throw new PrismLangException("Duplicated identifier \"" + s + "\"", constantList
						.getConstantNameIdent(i));
			} else {
				constantIdents.add(s);
			}
		}
	}

	// check variable names

	private void checkVarNames() throws PrismLangException
	{
		int i, j, n, m;
		Module module;
		String s;

		// compile list of all var names
		// and check as we go through

		// globals
		n = getNumGlobals();
		for (i = 0; i < n; i++) {
			s = getGlobal(i).getName();
			if (isIdentUsed(s)) {
				throw new PrismLangException("Duplicated identifier \"" + s + "\"", getGlobal(i));
			} else {
				varIdents.add(s);
				varDecls.add(getGlobal(i));
				varNames.add(s);
				varTypes.add(getGlobal(i).getType());
			}
		}

		// locals
		n = modules.size();
		for (i = 0; i < n; i++) {
			module = getModule(i);
			m = module.getNumDeclarations();
			for (j = 0; j < m; j++) {
				s = module.getDeclaration(j).getName();
				if (isIdentUsed(s)) {
					throw new PrismLangException("Duplicated identifier \"" + s + "\"", module.getDeclaration(j));
				} else {
					varIdents.add(s);
					varDecls.add(module.getDeclaration(j));
					varNames.add(s);
					varTypes.add(module.getDeclaration(j).getType());
				}
			}
		}

		// check there is at least one variable
		if (varNames.size() == 0) {
			throw new PrismLangException("There must be at least one variable");
		}
	}

	// Check there are no duplicate names labelling reward structs

	private void checkRewardStructNames() throws PrismLangException
	{
		int i, n;
		String s;
		HashSet<String> names = new HashSet<String>();
		n = getNumRewardStructs();
		for (i = 0; i < n; i++) {
			s = getRewardStruct(i).getName();
			if (s != null && !"".equals(s))
				if (!names.add(s))
					throw new PrismLangException("Duplicated reward structure name \"" + s + "\"", getRewardStruct(i));
		}
	}

	/**
	 * Get  a list of constants in the model that are undefined
	 * ("const int x;" rather than "const int x = 1;") 
	 */
	public Vector<String> getUndefinedConstants()
	{
		return constantList.getUndefinedConstants();
	}

	/**
	 * Set values for *all* undefined constants and then evaluate all constants.
	 * If there are no undefined constants, {@code someValues} can be null.
	 * Undefined constants can be subsequently redefined to different values with the same method.
	 * The current constant values (if set) are available via {@link #getConstantValues()}. 
	 * Calling this method also triggers some additional semantic checks
	 * that can only be done once constant values have been specified.
	 */
	public void setUndefinedConstants(Values someValues) throws PrismLangException
	{
		constantValues = constantList.evaluateConstants(someValues, null);
                for(Object module : modules) {
                    ((Module) module).setUndefinedConstants(someValues, constantValues);
                }
		semanticCheckAfterConstants(this, null);
	}

	/**
	 * Check if {@code name} is a *defined* constant in the model,
	 * i.e. a constant whose value was *not* left unspecified.
	 */
	public boolean isDefinedConstant(String name)
	{
		return constantList.isDefinedConstant(name);
	}
	
	/**
	 * Get access to the values for all constants in the model, including the 
	 * undefined constants set previously via the method {@link #setUndefinedConstants(Values)}.
	 * Until they are set for the first time, this method returns null.  
	 */
	public Values getConstantValues()
	{
		return constantValues;
	}

	/**
	 * Create a State object representing the default initial state of this model.
	 * If there are potentially multiple initial states (because the model has an 
	 * init...endinit specification), this method returns null;
	 * Assumes that values for constants have been provided for the model.
	 * Note: This method replaces the old getInitialValues() method,
	 * since State objects are now preferred to Values objects for efficiency.
	 */
	public State getDefaultInitialState() throws PrismLangException
	{
		int i, j, count, n, n2;
		Module module;
		Declaration decl;
		State initialState;
		Object initialValue;

		if (initStates != null) {
			return null;
		}

		// Create State object
		initialState = new State(getNumVars());
		// Then add values for all globals and all locals, in that order
		count = 0;
		n = getNumGlobals();
		for (i = 0; i < n; i++) {
			decl = getGlobal(i);
			initialValue = decl.getStartOrDefault().evaluate(constantValues);
			initialValue = getGlobal(i).getType().castValueTo(initialValue);
			initialState.setValue(count++, initialValue);
		}
		n = getNumModules();
		for (i = 0; i < n; i++) {
			module = getModule(i);
			n2 = module.getNumDeclarations();
			for (j = 0; j < n2; j++) {
				decl = module.getDeclaration(j);
				initialValue = decl.getStartOrDefault().evaluate(constantValues);
				initialValue = module.getDeclaration(j).getType().castValueTo(initialValue);
				initialState.setValue(count++, initialValue);
			}
		}

		return initialState;
	}

	/**
	 * Create a Values object representing the default initial state of this model.
	 * Deprecated: Use getDefaultInitialState() instead
	 * (or new Values(getDefaultInitialState(), modulesFile)).
	 */
	public Values getInitialValues() throws PrismLangException
	{
		int i, j, n, n2;
		Module module;
		Declaration decl;
		Values values;
		Object initialValue;

		if (initStates != null) {
			throw new PrismLangException("There are multiple initial states");
		}

		// set up variable list
		values = new Values();

		// first add all globals
		n = getNumGlobals();
		for (i = 0; i < n; i++) {
			decl = getGlobal(i);
			initialValue = decl.getStartOrDefault().evaluate(constantValues);
			initialValue = getGlobal(i).getType().castValueTo(initialValue);
			values.addValue(decl.getName(), initialValue);
		}
		// then add all module variables
		n = getNumModules();
		for (i = 0; i < n; i++) {
			module = getModule(i);
			n2 = module.getNumDeclarations();
			for (j = 0; j < n2; j++) {
				decl = module.getDeclaration(j);
				initialValue = decl.getStartOrDefault().evaluate(constantValues);
				initialValue = module.getDeclaration(j).getType().castValueTo(initialValue);
				values.addValue(decl.getName(), initialValue);
			}
		}

		return values;
	}

	/**
	 * Recompute all information about variables.
	 * More precisely... TODO
	 * Note: This does not re-compute the list of all identifiers used. 
	 */
	public void recomputeVariableinformation() throws PrismLangException
	{
		int i, n;

		// Recompute lists of all variables and types
		varDecls = new Vector<Declaration>();
		varNames = new Vector<String>();
		varTypes = new Vector<Type>();
		// Globals
		for (Declaration decl : globals) {
			varDecls.add(decl);
			varNames.add(decl.getName());
			varTypes.add(decl.getType());
		}
		// Locals
		n = modules.size();
		for (i = 0; i < n; i++) {
			for (Declaration decl : getModule(i).getDeclarations()) {
				varDecls.add(decl);
				varNames.add(decl.getName());
				varTypes.add(decl.getType());
			}
		}
		// Find all instances of variables, replace identifiers with variables.
		// Also check variables valid, store indices, etc.
		findAllVars(varNames, varTypes);
	}

	/**
	 * Create a VarList object storing information about all variables in this model.
	 * Assumes that values for constants have been provided for the model.
	 * Also performs various syntactic checks on the variables.   
	 */
	public VarList createVarList() throws PrismLangException
	{
		return new VarList(this);
	}

	// Methods required for ASTElement:

	/**
	 * Visitor method.
	 */
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}

	/**
	 * Convert to string.
	 */
	public String toString()
	{
		String s = "", tmp;
		int i, n;

		s += modelType.toString().toLowerCase() + "\n\n";

		tmp = "" + formulaList;
		if (tmp.length() > 0)
			tmp += "\n";
		s += tmp;

		tmp = "" + labelList;
		if (tmp.length() > 0)
			tmp += "\n";
		s += tmp;

		tmp = "" + constantList;
		if (tmp.length() > 0)
			tmp += "\n";
		s += tmp;

		n = getNumGlobals();
		for (i = 0; i < n; i++) {
			s += "global " + getGlobal(i) + ";\n";
		}
		if (n > 0) {
			s += "\n";
		}

		for (i = 0; i < modules.size() - 1; i++) {
			s += modules.elementAt(i) + "\n\n";
		}
		s += modules.elementAt(modules.size() - 1) + "\n";

		if (systemDefn != null) {
			s += "\nsystem " + systemDefn + " endsystem\n";
		}

		n = getNumRewardStructs();
		for (i = 0; i < n; i++) {
			s += "\n" + getRewardStruct(i);
		}

		if (initStates != null) {
			s += "\n" + "init " + initStates + " endinit\n";
		}

		return s;
	}

	/**
	 * Perform a deep copy.
	 */
	@SuppressWarnings("unchecked")
	public ASTElement deepCopy()
	{
		int i, n;
		ModulesFile ret = new ModulesFile();
		
		// Copy ASTElement stuff
		ret.setPosition(this);
		// Copy type
		ret.setModelType(modelType);
		// Deep copy main components
		ret.setFormulaList((FormulaList) formulaList.deepCopy());
		ret.setLabelList((LabelList) labelList.deepCopy());
		ret.setConstantList((ConstantList) constantList.deepCopy());
		n = getNumGlobals();
		for (i = 0; i < n; i++) {
			ret.addGlobal((Declaration) getGlobal(i).deepCopy());
		}
		n = getNumModules();
		for (i = 0; i < n; i++) {
			ret.addModule((Module) getModule(i).deepCopy());
		}
		if (systemDefn != null)
			ret.setSystemDefn(systemDefn.deepCopy());
		n = getNumRewardStructs();
		for (i = 0; i < n; i++) {
			ret.addRewardStruct((RewardStruct) getRewardStruct(i).deepCopy());
		}
		if (initStates != null)
			ret.setInitialStates(initStates.deepCopy());
		// Copy other (generated) info
		ret.formulaIdents = (formulaIdents == null) ? null : (Vector<String>)formulaIdents.clone();
		ret.constantIdents = (constantIdents == null) ? null : (Vector<String>)constantIdents.clone();
		ret.varIdents = (varIdents == null) ? null : (Vector<String>)varIdents.clone();
		ret.moduleNames = (moduleNames == null) ? null : moduleNames.clone();
		ret.synchs = (synchs == null) ? null : (Vector<String>)synchs.clone();
		ret.synchsOfFDToModule = (synchsOfFDToModule == null) ? null : new HashMap<String,Integer>(synchsOfFDToModule);
		ret.synchsOfFDToFD = (synchsOfFDToFD == null) ? null : new HashMap<String,Integer>(synchsOfFDToFD);
		if (varDecls != null) {
			ret.varDecls = new Vector<Declaration>();
			for (Declaration d : varDecls)
				ret.varDecls.add((Declaration) d.deepCopy());
		}
		ret.varNames = (varNames == null) ? null : (Vector<String>)varNames.clone();
		ret.varTypes = (varTypes == null) ? null : (Vector<Type>)varTypes.clone();
		ret.constantValues = (constantValues == null) ? null : new Values(constantValues);
		
		return ret;
	}
}

// ------------------------------------------------------------------------------
